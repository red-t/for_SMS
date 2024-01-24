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
    ('startIdx',        np.int32),
    ('endIdx',          np.int32),
    ('teTid',           np.int32),
])

TeAlignmentDt = np.dtype([
    ('segIdx',      np.int32),
    ('AlnScore',    np.int32),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
    ('mapLen',      np.int32),
    ('divergence',  np.float32),
    ('flag',        np.int16),
    ('teTid',       np.int32),
])

ClusterDt = np.dtype([
    ('tid',                 np.int32),
    ('refStart',            np.int32),
    ('refEnd',              np.int32),
    ('idx',                 np.int32),
    ('startIdx',            np.int32),
    ('endIdx',              np.int32),
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
    ('insEnd',              np.int32),
    ('insStart',            np.int32),
    ('flag',                np.uint16),
])


##########################
### Construct SegArray ###
##########################
cdef object getSegArray(BamFile genomeBam, Args args):
    cdef int returnValue, numSeg=0, maxNum=9900
    cdef object segArray = np.zeros(10000, dtype=SegmentDt)
    cdef Segment[::1] segArrayView = segArray
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    
    while True:
        returnValue = iterator.cnext1()
        if returnValue < 0:
            del iterator
            segArray.resize((numSeg,), refcheck=False)
            return segArray

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


cdef updateSegArrayByTe(Segment[::1] segArray, Args args):
    cdef BamFile teBam = BamFile("tmp_build/all_seg_{}.bam".format(args.tid), "rb", args.numThread)
    cdef Iterator iterator = Iterator(teBam)
    cdef object teArray = getTeArray(iterator)
    cdef TeAlignment[::1] teArrayView = teArray
    cdef int numTeTid = teBam.header.n_targets

    teArray.sort(order=['segIdx', 'queryStart'])
    
    cdef object teTidCountTable = np.zeros(numTeTid, dtype=np.int32)
    cdef int[::1] teTidCountTableView = teTidCountTable
    cdef int i

    for i in range(teArrayView.shape[0]):
        updateSegByTeArray(&segArray[0], &teArrayView[0], i)

    for i in range(segArray.shape[0]):
        if segArray[i].numTeAlignment:
            countTeTids(&segArray[i], &teArrayView[0], &teTidCountTableView[0], numTeTid)
            segArray[i].teTid = np.argmax(teTidCountTable)

    del iterator; teBam.close(); del teBam; del teTidCountTable


#########################
### Construct TeArray ###
#########################
cdef mapSegToTE(str teFn, Args args):
    cdef int exitCode
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t {} -aY {} tmp_build/all_seg_{}.fa | " \
                   "samtools view -@ {} -bhS -o tmp_build/all_seg_{}.bam -".format(args.numThread, teFn, args.tid, args.numThread, args.tid)

    process = Popen([cmd], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitCode = process.wait()
    if exitCode != 0:
        raise Exception("Error: minimap2 failed for tmp_build/all_seg_{}.fa".format(args.tid))


cdef object getTeArray(Iterator iterator):
    cdef int returnValue, numTeAln=0, maxNum=9900
    cdef object teArray = np.zeros(10000, dtype=TeAlignmentDt)
    cdef TeAlignment[::1] teArrayView = teArray
    
    while True:
        returnValue = iterator.cnext2()
        if returnValue < 0:
            teArray.resize((numTeAln,), refcheck=False)
            return teArray

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numTeAln >= maxNum:
            maxNum = teArray.shape[0] + 10000
            teArray.resize((maxNum,), refcheck=False)
            teArrayView = teArray
            maxNum -= 100

        fillTeArray(iterator.bamRcord, &teArrayView[numTeAln])
        numTeAln += 1


##########################
### Construct CltArray ###
##########################
cdef object getCltArray(Segment[::1] segArray, Args args):
    cdef int startIdx=0, endIdx, numClt=0, maxNum=9900
    cdef object cltArray = np.zeros(10000, dtype=ClusterDt)
    cdef Cluster[::1] cltArrayView = cltArray
    
    while startIdx < segArray.shape[0]:
        if overhangIsShort(&segArray[startIdx], args.minOverhang):
            startIdx += 1; continue

        if numClt > maxNum:
            maxNum = cltArray.shape[0] + 10000
            cltArray.resize((maxNum,), refcheck=False)
            cltArrayView = cltArray
            maxNum -= 100
        
        # Initialize numClt-th cluster
        cltArrayView[numClt].tid = args.tid
        cltArrayView[numClt].refStart = segArray[startIdx].refPosition - 1
        cltArrayView[numClt].refEnd = segArray[startIdx].refPosition + args.maxDistance
        cltArrayView[numClt].idx = numClt
        cltArrayView[numClt].startIdx = startIdx

        endIdx = startIdx + 1
        while endIdx < segArray.shape[0]:
            if overhangIsShort(&segArray[endIdx], args.minOverhang):
                endIdx += 1; continue
            if segArray[endIdx].refPosition > cltArrayView[numClt].refEnd:
                break

            cltArrayView[numClt].refEnd = segArray[endIdx].refPosition + args.maxDistance
            endIdx += 1
        
        cltArrayView[numClt].endIdx = endIdx
        cltArrayView[numClt].refEnd = cltArrayView[numClt].refEnd - args.maxDistance
        startIdx = endIdx; numClt += 1
    
    cltArray.resize((numClt,), refcheck=False)
    return cltArray


cdef updateCltArray(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBam, Args args):
    
    cdef BamFile teBam = BamFile("tmp_build/all_seg_{}.bam".format(args.tid), "rb", 1)
    cdef object teTidCountTable = np.zeros(teBam.header.n_targets, dtype=np.int32)
    cdef int[::1] teTidCountTableView = teTidCountTable
    cdef int i
    
    args.numTeTid = teBam.header.n_targets
    args.teTidCountTable = &teTidCountTableView[0]
    args.genomeBam = genomeBam.htsFile
    args.firstBamRecord = bam_init1()
    args.secondBamRecord = bam_init1()

    teBam.close(); del teBam
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
cdef outPut(object cltArray, Segment[::1] segArray, BamFile genomeBam, Args args):

    cdef bytes qnameBytes
    cdef bytes chromBytes = sam_hdr_tid2name(genomeBam.header, args.tid)
    cdef str cltId, chrom = chromBytes.decode()
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
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
        for j in range(cltArray[i]['startIdx'], cltArray[i]['endIdx']):
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
    

#####################
### Build Cluster ###
#####################
cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid):
    # 1. Construct segments
    cdef Args args = newArgs(tid, bgDiv, bgDepth, bgReadLen, cmdArgs)
    cdef BamFile genomeBam = BamFile(cmdArgs.genomeBamFilePath, "rb", cmdArgs.numThread)
    cdef object segArray = getSegArray(genomeBam, args)

    # 2. Compute segment features
    cdef const char *chrom = sam_hdr_tid2name(genomeBam.header, tid)
    args.repeatAiList = newAiList(cmdArgs.repeatFn, chrom)
    args.gapAiList = newAiList(cmdArgs.gapFn, chrom)

    updateSegArray(segArray, args)
    segArray.sort(order='refPosition')
    ouputAllSegSeqs(segArray, genomeBam, args)

    mapSegToTE(cmdArgs.teFn, args)
    updateSegArrayByTe(segArray, args)
    
    # 3. Construct cluster
    cdef tidToCltData = {}
    cdef object cltArray = getCltArray(segArray, args)

    # 4. Compute cluster features
    updateCltArray(cltArray, segArray, genomeBam, args)

    # 5. Filter clusters
    args.blackAiList = newAiList(cmdArgs.blackListPath, chrom)
    filterByBlacklist(cltArray, args)
    cltArray = filterByModel(cltArray, cmdArgs)

    # 6. Output cluster seqs for assembling
    outputGermCltSeqs(cltArray, segArray, genomeBam, args)

    tidToCltData[tid] = (cltArray, segArray)
    genomeBam.close(); del genomeBam
    return tidToCltData


##############################
### Get High-Qual Clusters ###
##############################
cdef object getHighQualClts(dict allCltData):
    cdef int i, tid, numClt = 0, maxNum = 1900
    cdef object highQualArray = np.zeros(2000, dtype=ClusterDt)
    cdef Cluster[::1] arrayView = highQualArray
    cdef Cluster[::1] cltArray

    for tid in allCltData.keys():
        cltArray = allCltData[tid][0]
        if cltArray.shape[0] == 0:
            continue
        
        for i in range(cltArray.shape[0]):
            if isLowQualClt(&cltArray[i]):
                continue
            
            if numClt >= maxNum:
                maxNum = highQualArray.shape[0] + 2000
                highQualArray.resize((maxNum,), refcheck=False)
                arrayView = highQualArray
                maxNum -= 100
            
            arrayView[numClt] = cltArray[i]
            numClt += 1

    highQualArray.resize((numClt,), refcheck=False)
    return highQualArray