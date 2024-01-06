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

##############################
### Construction Functions ###
##############################
cdef AiList* newAiList(str filePath, const char *chrom):
    cdef bytes filePathBytes = filePath.encode()
    cdef AiList *aiList = initAiList()

    readBED(aiList, filePathBytes, chrom)
    constructAiList(aiList, 20)
    return aiList


cdef Args newArgs(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs):
    cdef Args args

    args.tid = tid
    args.bgDiv = bgDiv
    args.bgDepth = bgDepth
    args.bgReadLen = bgReadLen
    args.numThread = cmdArgs.numThread
    args.minSegLen = cmdArgs.minSegLen
    args.maxDistance = cmdArgs.maxDistance
    args.minOverhang = cmdArgs.minOverhang
    return args


##########################
### Construct SegArray ###
##########################
cdef object getSegArray(BamFile genomeBamFile, Args args):

    cdef object segArray = np.zeros(10000, dtype=SegmentDt)
    cdef object template = np.zeros(10000, dtype=SegmentDt)
    cdef Segment[::1] segArrayView = segArray
    cdef int maxNumSeg = segArrayView.shape[0] - 20
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    # cdef BamFile outputBamFile = BamFile("tmp_candidates_alignments.{}.bam".format(args.tid), "wb", 5, genomeBamFile)
    cdef int returnValue, numSeg=0
    
    while True:
        returnValue = iterator.cnext1()
        if returnValue < 0:
            # outputBamFile.close(); del outputBamFile
            del template; del iterator
            return segArray[:numSeg]

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numSeg > maxNumSeg:
            segArray = np.concatenate((segArray, template))
            segArrayView = segArray
            maxNumSeg = segArrayView.shape[0] - 20

        returnValue = fillSegmentArray(iterator.bamRcord, &segArrayView[numSeg], iterator.offset, args.minSegLen)
        numSeg += returnValue
        # if returnValue > 0:
            # outputBamFile.write(iterator.bamRcord)


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


##########################
### Construct TeArray ###
##########################
cdef ouputSegmentSeqs(Segment[::1] segArray, BamFile genomeBamFile, Args args):

    cdef str outputFileName = "tmp.all_supp_reads.{}.fa".format(args.tid)
    cdef BamFile outputFasta = BamFile(outputFileName, "wF", args.numThread, genomeBamFile)
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, returnValue

    for i in range(segArray.shape[0]):
        returnValue = iterator.cnext3(segArray[i].fileOffset)
        trimSegment(iterator.bamRcord, destRecord, i, segArray[i].queryStart, segArray[i].queryEnd)
        outputFasta.write(destRecord)
    
    bam_destroy1(destRecord); outputFasta.close(); del outputFasta; del iterator


cdef mapByMinimap2(str reference, Args args):

    cdef int retval
    cdef str command = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t {} -aY {} tmp.all_supp_reads.{}.fa | " \
                      "samtools view -@ {} -bhS -o tmp.all_supp_reads.{}.bam -".format(args.numThread, reference, args.tid, args.numThread, args.tid)

    process = Popen([command], stderr=DEVNULL, shell=True, executable='/bin/bash')
    retval = process.wait()
    if retval != 0:
        raise Exception("Error: minimap2 failed for tmp.all_supp_reads.{}.fa".format(args.tid))


cdef object getTeArray(Iterator iterator):

    cdef object teArray = np.zeros(10000, dtype=TeAlignmentDt)
    cdef object template  = np.zeros(10000, dtype=TeAlignmentDt)
    cdef TeAlignment[::1] teArrayView = teArray
    cdef int maxTeAlignments = teArrayView.shape[0] - 20
    cdef int returnValue, numTeAlignments=0
    
    while True:
        returnValue = iterator.cnext2()
        if returnValue < 0:
            del template
            return teArray[:numTeAlignments]

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numTeAlignments > maxTeAlignments:
            teArray = np.concatenate((teArray, template))
            teArrayView = teArray
            maxTeAlignments = teArrayView.shape[0] - 20

        fillTeArray(iterator.bamRcord, &teArrayView[numTeAlignments])
        numTeAlignments += 1


##########################
### Construct CltArray ###
##########################
cdef object getCltArray(Segment[::1] segArray, Args args):

    cdef object cltArray = np.zeros(10000, dtype=ClusterDt)
    cdef object template = np.zeros(10000, dtype=ClusterDt)
    cdef Cluster[::1] cltArrayView = cltArray
    cdef int maxNumClt = cltArrayView.shape[0] - 20
    cdef int numClt=0, start=0, end
    
    while start < segArray.shape[0]:
        if overhangIsShort(&segArray[start], args.minOverhang):
            start += 1; continue

        if numClt > maxNumClt:
            cltArray = np.concatenate((cltArray, template))
            cltArrayView = cltArray
            maxNumClt = cltArrayView.shape[0] - 20
        
        # initialize numClt-th cluster
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
    
    del template; return cltArray[:numClt]


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


######################
### Local Assembly ###
######################
cdef assembleClusters(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args):
    if cltArray.shape[0] == 0:
        return

    outputGermSeqs(cltArray, segArray, genomeBamFile, args)
    # wtdbg2Assemble(cltArray)
    outputSomaSeqs(cltArray, segArray, genomeBamFile, args)

cdef outputGermSeqs(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args):
    cdef str outputFileName
    cdef BamFile outputFasta
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, j

    for i in range(cltArray.shape[0]):
        if isLowQualClt(&cltArray[i]) or isSomaClt(&cltArray[i]):
            continue

        outputFileName = "tmp.{}_{}.fa".format(args.tid, i)
        outputFasta = BamFile(outputFileName, "wF", args.numThread, genomeBamFile)
        for j in range(cltArray[i].startIndex, cltArray[i].endIndex):
            if overhangIsShort(&segArray[j], args.minOverhang):
                continue
            outputSingleSeq(segArray, outputFasta, iterator, destRecord, j)
        
        outputFasta.close()

    bam_destroy1(destRecord); del iterator

cdef outputSingleSeq(Segment[::1]segArray, BamFile outputFasta, Iterator iterator, bam1_t *destRecord, int j, int flankSize=3000):
    cdef int start, end, returnValue

    returnValue = iterator.cnext3(segArray[j].fileOffset)
    getTrimRegion(&segArray[j], &start, &end, flankSize)
    trimSegment(iterator.bamRcord, destRecord, j, start, end)
    outputFasta.write(destRecord)

cdef wtdbg2Assemble(Cluster[::1] cltArray):
    pass

cdef outputSomaSeqs(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args):
    cdef str outputFileName
    cdef BamFile outputFasta
    cdef Iterator iterator = Iterator(genomeBamFile, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, j

    for i in range(cltArray.shape[0]):
        if isLowQualClt(&cltArray[i]):
            continue
        
        # Skip successfully assembled clusters
        outputFileName = "tmp.{}_{}_assembled.fa".format(args.tid, i)
        if os.path.isfile(outputFileName):
            continue
        
        outputFasta = BamFile(outputFileName, "wF", args.numThread, genomeBamFile)
        j = getOuputSegIndex(&cltArray[i], &segArray[0], args)
        outputSingleSeq(segArray, outputFasta, iterator, destRecord, j)
        outputFasta.close()

    bam_destroy1(destRecord); del iterator
    

############
### Main ###
############
cpdef dict buildCluster(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs):
    # 1. construct segments
    cdef Args args = newArgs(tid, bgDiv, bgDepth, bgReadLen, cmdArgs)
    cdef BamFile genomeBamFile = BamFile(cmdArgs.genomeBamFilePath, "rb", cmdArgs.numThread)
    cdef object segArray = getSegArray(genomeBamFile, args)
    cdef Segment[::1] segArrayView = segArray

    # 2. compute segment features
    cdef const char *chrom = sam_hdr_tid2name(genomeBamFile.header, tid)
    args.repeatAiList = newAiList(cmdArgs.repeatPath, chrom)
    args.gapAiList = newAiList(cmdArgs.gapPath, chrom)

    updateSegArray(segArrayView, args)
    segArray.sort(order='refPosition')
    ouputSegmentSeqs(segArrayView, genomeBamFile, args)

    mapByMinimap2(cmdArgs.referenceTe, args)
    teArray = updateSegArrayByTe(segArrayView, args)
    
    # 3. construct cluster
    cdef ResultDict = {}
    cdef object cltArray = getCltArray(segArray, args)
    cdef Cluster[::1] cltArrayView = cltArray

    # 4. compute cluster features
    updateCltArray(cltArrayView, segArrayView, genomeBamFile, args)

    # 5. filter cluster
    args.blackAiList = newAiList(cmdArgs.blackListPath, chrom)
    filterByBlacklist(cltArrayView, args)
    cltArray = filterByModel(cltArray, cmdArgs)

    # 6. local assembly
    cltArrayView = cltArray
    assembleClusters(cltArrayView, segArrayView, genomeBamFile, args)

    # 6. output
    outPut(cltArray, segArrayView, genomeBamFile, args)
    genomeBamFile.close(); del genomeBamFile

    ResultDict[tid] = (cltArray, segArray, teArray)
    return ResultDict