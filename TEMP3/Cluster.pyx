import os
import numpy as np
import pandas as pd
from subprocess import Popen, DEVNULL
from autogluon.tabular import TabularPredictor

##################
### Data Types ###
##################
SegmentDt = np.dtype([
    ('flag',            np.uint16),
    ('mapQual',         np.uint8),
    ('queryStart',      np.int32),
    ('queryEnd',        np.int32),
    ('refPosition',     np.int32),
    ('segType',         np.uint8),
    ('alnType',         np.uint8),
    ('fileOffset',      np.int64),
    ('alnRefStart',     np.int32),
    ('alnRefEnd',       np.int32),
    ('order',           np.uint8),
    ('numSeg',          np.uint8),
    ('overhang',        np.int32),
    ('matchLen',        np.int32),
    ('readLen',         np.int32),
    ('alnLocationType', np.uint8),
    ('numTeAlignment',  np.uint8),
    ('sumQueryMapLen',  np.int32),
    ('sumAlnScore',     np.float32),
    ('sumDivergence',   np.float32),
    ('directionFlag',   np.uint16),
    ('startIndex',      np.int32),
    ('endIndex',        np.int32),
    ('teTid',           np.int32),
])

TeAlignmentDt = np.dtype([
    ('segIndex',    np.int32),
    ('AlnScore',    np.int32),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
    ('mapLen',      np.int32),
    ('divergence',  np.float32),
    ('flag',        np.int16),
    ('teTid',       np.int32),
])

ClusterDt = np.dtype([
    ('refStart',            np.int32),
    ('refEnd',              np.int32),
    ('startIndex',          np.int32),
    ('endIndex',            np.int32),
    ('numSeg',              np.float32),
    ('directionFlag',       np.uint16),
    ('cltType',             np.uint8),
    ('locationType',        np.uint8),
    ('numSegType',          np.uint8),
    ('entropy',             np.float32),
    ('balanceRatio',        np.float32),
    ('lowMapQualFrac',      np.float32),
    ('dualClipFrac',        np.float32),
    ('alnFrac1',            np.float32),
    ('alnFrac2',            np.float32),
    ('alnFrac4',            np.float32),
    ('alnFrac8',            np.float32),
    ('alnFrac16',           np.float32),
    ('meanMapQual',         np.float32),
    ('meanAlnScore',        np.float32),
    ('meanQueryMapFrac',    np.float32),
    ('meanDivergence',      np.float32),
    ('bgDiv',               np.float32),
    ('bgDepth',             np.float32),
    ('bgReadLen',           np.float32),
    ('teAlignedFrac',       np.float32),
    ('teTid',               np.int32),
    ('isInBlacklist',       np.uint8),
    ('probability',         np.float32),
])

#######################
### Construc AiList ###
#######################
cdef AiList* newAiList(str filePath, const char *chrom):
    cdef bytes filePathBytes = filePath.encode()
    cdef AiList *aiList = initAiList()

    readBED(aiList, filePathBytes, chrom)
    constructAiList(aiList, 20)
    return aiList


##########################
### Construct SegArray ###
##########################
cdef object getSegArray(BamFile genomeBamFile, Args args):
    cdef int returnValue, numSeg=0, maxNum=9900
    cdef object segArray = np.zeros(10000, dtype=SegmentDt)
    cdef Segment[::1] segArrayView = segArray
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    
    while True:
        returnValue = iterator.cnext1()
        if returnValue < 0:
            del iterator
            return segArray[:numSeg]

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numSeg >= maxNum:
            maxNum = segArrayView.shape[0] + 10000
            segArray.resize((maxNum,), refcheck=False)
            segArrayView = segArray
            maxNum -= 100

        returnValue = fillSegmentArray(iterator.bamRcord, &segArrayView[numSeg], iterator.offset, args.minSegLen)
        numSeg += returnValue


cdef updateSegArray(Segment[::1] segArray, Args args):
    cdef int i
    for i in range(segArray.shape[0]):
        updateSegment(&segArray[i], args.repeatAiList, args.gapAiList)


cdef object updateSegArrayByTe(Segment[::1] segArray, Args args):
    cdef BamFile teBamFile = BamFile("tmp.all_supp_reads.{}.bam".format(args.tid), "rb", args.numThread)
    cdef Iterator iterator = Iterator(teBamFile)
    cdef object teArray = getTeArray(iterator)
    cdef TeAlignment[::1] teArrayView = teArray
    cdef int numTeTid = teBamFile.header.n_targets

    teArray.sort(order=['segIndex', 'queryStart'])
    
    cdef object teTidCountTable = np.zeros(numTeTid, dtype=np.int32)
    cdef int[::1] teTidCountTableView = teTidCountTable
    cdef int i

    for i in range(teArrayView.shape[0]):
        updateSegByTeArray(&segArray[0], &teArrayView[0], i)

    for i in range(segArray.shape[0]):
        if segArray[i].numTeAlignment:
            countTeTids(&segArray[i], &teArrayView[0], &teTidCountTableView[0], numTeTid)
            segArray[i].teTid = np.argmax(teTidCountTable)

    del iterator; teBamFile.close(); del teBamFile; del teTidCountTable
    return teArray


#########################
### Construct TeArray ###
#########################
cdef mapByMinimap2(str reference, Args args):
    cdef int returnValue
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t {} -aY {} tmp.all_supp_reads.{}.fa | " \
                   "samtools view -@ {} -bhS -o tmp.all_supp_reads.{}.bam -".format(args.numThread, reference, args.tid, args.numThread, args.tid)

    process = Popen([cmd], stderr=DEVNULL, shell=True, executable='/bin/bash')
    returnValue = process.wait()
    if returnValue != 0:
        raise Exception("Error: minimap2 failed for tmp.all_supp_reads.{}.fa".format(args.tid))


cdef object getTeArray(Iterator iterator):
    cdef int returnValue, numTeAlignments=0, maxNum=9900
    cdef object teArray = np.zeros(10000, dtype=TeAlignmentDt)
    cdef TeAlignment[::1] teArrayView = teArray
    
    while True:
        returnValue = iterator.cnext2()
        if returnValue < 0:
            return teArray[:numTeAlignments]

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numTeAlignments >= maxNum:
            maxNum = teArray.shape[0] + 10000
            teArray.resize((maxNum,), refcheck=False)
            teArrayView = teArray
            maxNum -= 100

        fillTeArray(iterator.bamRcord, &teArrayView[numTeAlignments])
        numTeAlignments += 1


##########################
### Construct CltArray ###
##########################
cdef object getCltArray(Segment[::1] segArray, Args args):
    cdef int start=0, end, numClt=0, maxNum=9900
    cdef object cltArray = np.zeros(10000, dtype=ClusterDt)
    cdef Cluster[::1] cltArrayView = cltArray
    
    while start < segArray.shape[0]:
        if overhangIsShort(&segArray[start], args.minOverhang):
            start += 1; continue

        if numClt > maxNum:
            maxNum = cltArray.shape[0] + 10000
            cltArray.resize((maxNum,), refcheck=False)
            cltArrayView = cltArray
            maxNum -= 100
        
        # Initialize numClt-th cluster
        cltArrayView[numClt].startIndex = start
        cltArrayView[numClt].refStart = segArray[start].refPosition - 1
        cltArrayView[numClt].refEnd = segArray[start].refPosition + args.maxDistance

        end = start + 1
        while end < segArray.shape[0]:
            if overhangIsShort(&segArray[end], args.minOverhang):
                end += 1; continue
            if segArray[end].refPosition > cltArrayView[numClt].refEnd:
                break

            cltArrayView[numClt].refEnd = segArray[end].refPosition + args.maxDistance
            end += 1
        
        cltArrayView[numClt].endIndex = end
        cltArrayView[numClt].refEnd = cltArrayView[numClt].refEnd - args.maxDistance
        start = end; numClt += 1
    
    return cltArray[:numClt]


cdef updateCltArray(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args):
    
    cdef BamFile teBamFile = BamFile("tmp.all_supp_reads.{}.bam".format(args.tid), "rb", 1)
    cdef object teTidCountTable = np.zeros(teBamFile.header.n_targets, dtype=np.int32)
    cdef int[::1] teTidCountTableView = teTidCountTable
    cdef int i
    
    args.numTeTid = teBamFile.header.n_targets
    args.teTidCountTable = &teTidCountTableView[0]
    args.genomeBamFile = genomeBamFile.htsFile
    args.firstBamRecord = bam_init1()
    args.secondBamRecord = bam_init1()

    teBamFile.close(); del teBamFile
    for i in range(cltArray.shape[0]):
        updateCluster(&cltArray[i], &segArray[0], args)

        if not isValidCandidate(&cltArray[i]):
            continue
        cltArray[i].teTid = np.argmax(teTidCountTable)
    
    bam_destroy1(args.firstBamRecord); bam_destroy1(args.secondBamRecord)
    destroyAiList(args.repeatAiList); destroyAiList(args.gapAiList); del teTidCountTable


######################
### Filter Cluster ###
######################
cdef filterByBlacklist(Cluster[::1] cltArray, Args args):
    cdef int i

    for i in range(cltArray.shape[0]):
        intersectBlackList(&cltArray[i], args)
    
    destroyAiList(args.blackAiList)


cdef object filterByModel(object cltArray, object cmdArgs):
    cdef object cltDf = pd.DataFrame(cltArray)
    
    filterGermByModel(cltDf, cmdArgs.germModelPath)
    filterSomaByModel(cltDf, cmdArgs.somaModelPath)

    return cltDf.to_records(index=False)

cdef filterGermByModel(object cltDf, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object germDf = cltDf.loc[(cltDf['cltType']==0) & (cltDf['teAlignedFrac']>=0.8) & (cltDf['isInBlacklist']==0)]

    probability = predictor.predict_proba(germDf)
    probability.columns = ['0', 'probability']
    cltDf.update(probability)

cdef filterSomaByModel(object cltDf, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object somaDf = cltDf.loc[(cltDf['cltType']>0) & (cltDf['teAlignedFrac']>=0.8) & (cltDf['isInBlacklist']==0)]

    probability = predictor.predict_proba(somaDf)
    probability.columns = ['0', 'probability']
    cltDf.update(probability)


##############
### Output ###
##############
cdef outPut(object cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args):

    cdef bytes qnameBytes
    cdef bytes chromBytes = sam_hdr_tid2name(genomeBamFile.header, args.tid)
    cdef str cltId, chrom = chromBytes.decode()
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    cdef list cltList
    cdef int i, j

    cltOutput = open('tmp_clt_{}.txt'.format(args.tid), 'w')
    segOutput = open('tmp_seg_{}.txt'.format(args.tid), 'w')

    for i in range(cltArray.shape[0]):
        ### output clt ###
        cltList = list(cltArray[i])

        # direction
        if cltList[5] == 1:
            cltList.insert(2, '+')
        elif cltList[5] == 2:
            cltList.insert(2, '-')
        else:
            cltList.insert(2, '*')
        
        # normalized numSeg
        cltList.insert(2, cltList[5])

        # cluster id
        cltId = str(args.tid) + '-' + str(i)
        cltList.insert(2, cltId)
        
        # chromosome
        cltList.insert(0, chrom)

        # write out clt
        cltList = [str(x) for x in cltList]
        cltOutput.write('\t'.join(cltList) + '\n')

        ### output seg ###
        for j in range(cltArray[i]['startIndex'], cltArray[i]['endIndex']):
            if segArray[j].overhang < args.minOverhang:
                continue

            # chromosome & cluster id
            cltList = [chrom, cltId]

            # refst
            cltList.append(str(segArray[j].alnRefStart))

            # qname
            iterator.cnext3(segArray[j].fileOffset)

            qnameBytes = bam_get_qname(iterator.bamRcord)
            cltList.append(qnameBytes.decode())

            # write out seg
            segOutput.write('\t'.join(cltList) + '\n')
    
    cltOutput.close(); segOutput.close(); del iterator
    

############
### Main ###
############
cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid):
    # 1. construct segments
    cdef Args args = newArgs(tid, bgDiv, bgDepth, bgReadLen, cmdArgs)
    cdef BamFile genomeBamFile = BamFile(cmdArgs.genomeBamFilePath, "rb", cmdArgs.numThread)
    cdef object segArray = getSegArray(genomeBamFile, args)

    # 2. compute segment features
    cdef const char *chrom = sam_hdr_tid2name(genomeBamFile.header, tid)
    args.repeatAiList = newAiList(cmdArgs.repeatPath, chrom)
    args.gapAiList = newAiList(cmdArgs.gapPath, chrom)

    updateSegArray(segArray, args)
    segArray.sort(order='refPosition')
    ouputAllSegSeqs(segArray, genomeBamFile, args)

    mapByMinimap2(cmdArgs.referenceTe, args)
    teArray = updateSegArrayByTe(segArray, args)
    
    # 3. construct cluster
    cdef chromCltData = {}
    cdef object cltArray = getCltArray(segArray, args)

    # 4. compute cluster features
    updateCltArray(cltArray, segArray, genomeBamFile, args)

    # 5. filter cluster
    args.blackAiList = newAiList(cmdArgs.blackListPath, chrom)
    filterByBlacklist(cltArray, args)
    cltArray = filterByModel(cltArray, cmdArgs)

    # 6. output sequences for local assembly
    outputGermCltSeqs(cltArray, segArray, genomeBamFile, args)

    chromCltData[tid] = (cltArray, segArray, teArray)
    genomeBamFile.close(); del genomeBamFile
    return chromCltData