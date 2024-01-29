import numpy as np
from .FileIO import outputSomaCltSeqs, outputRefFlank
from .Cluster import buildCluster
from .Assemble import assembleCluster
from .Annotate import annotateCluster
from concurrent.futures import ProcessPoolExecutor, as_completed


###########################
### Get Background Info ###
###########################
cdef dict getBackgroundInfo(str genomeBamFilePath, int numThread):
    cdef BamFile genomeBam = BamFile(genomeBamFilePath, "rb", numThread)
    cdef int i, tid=0, maxChromLen=0

    for i in range(genomeBam.header.n_targets):
        if maxChromLen < sam_hdr_tid2len(genomeBam.header, i):
            maxChromLen = sam_hdr_tid2len(genomeBam.header, i)
            tid = i
    
    cdef Iterator iterator = Iterator(genomeBam, tid)
    cdef int numAln = 0, maxNumAln = 499900
    cdef object readLenArray = np.zeros(500000, dtype=np.int32)
    cdef object divArray = np.zeros(500000, dtype=np.float32)
    cdef float[::1] divArrayView = divArray
    cdef int[::1] readLenArrayView = readLenArray
    cdef int retValue, alnLen
    cdef int64_t sumAlnLen = 0
    cdef float divergence
    cdef dict bgInfo = {}

    while True:
        retValue = iterator.cnext1()
        if retValue < 0:
            bgInfo["numChrom"] = genomeBam.header.n_targets
            bgInfo["bgDiv"] = np.mean(divArray[:numAln])
            bgInfo["bgDepth"] = float(sumAlnLen) / maxChromLen
            bgInfo["bgReadLen"] = np.median(readLenArray[:numAln])

            genomeBam.close(); del iterator; del genomeBam; return bgInfo

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numAln >= maxNumAln:
            maxNumAln = divArray.shape[0] + 500000
            divArray.resize((maxNumAln,), refcheck=False)
            readLenArray.resize((maxNumAln,), refcheck=False)
            divArrayView = divArray
            readLenArrayView = readLenArray
            maxNumAln -= 100
            
        getMapLenAndDiv(&alnLen, &divergence, iterator.bamRcord)
        sumAlnLen += alnLen
        divArrayView[numAln] = divergence
        readLenArrayView[numAln] = iterator.bamRcord.core.l_qseq
        numAln += 1


#######################
### Parallel Module ###
#######################
cdef tuple divideTask(int numTask, int poolSize):
    cdef int i, taskSize = int(numTask / poolSize)
    cdef list startList = []

    if taskSize == 0:
        taskSize = 1
    
    for i in range(0, numTask, taskSize):
        startList.append(i)
    
    return taskSize, startList


cpdef object runInParallel(object cmdArgs):    
    cdef set subProcTup
    cdef dict allCltData = {}, tidToCltData, bgInfo
    cdef object subProc, assembleArray, retValue
    cdef int startIdx, taskSize
    cdef list startList

    # 1. Get Background Info
    bgInfo = getBackgroundInfo(cmdArgs.genomeBamFilePath, cmdArgs.numThread)
    print("bg Divergence: {}\nbg coverage: {}\nbg readlen: {}" \
          "".format(bgInfo["bgDiv"], bgInfo["bgDepth"], bgInfo["bgReadLen"]))

    with ProcessPoolExecutor(max_workers=cmdArgs.numProcess) as executor:
        # 2. Build Cluster
        subProcTup = set([executor.submit(buildCluster, \
                                          bgInfo["bgDiv"], \
                                          bgInfo["bgDepth"], \
                                          bgInfo["bgReadLen"], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
        for subProc in as_completed(subProcTup):
            tidToCltData = subProc.result()
            allCltData = {**allCltData, **tidToCltData}
        
        # 3. Local Assembly
        highQualArray = getHighQualClts(allCltData)
        taskSize, startList = divideTask(highQualArray.shape[0], cmdArgs.numProcess)
        subProcTup = set([executor.submit(assembleCluster, \
                                          highQualArray, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs.numThread) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 4. Output sequence for clusters without assembly
        subProcTup = set([executor.submit(outputSomaCltSeqs, \
                                          allCltData[tid][0], \
                                          allCltData[tid][1], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 5. Output reference flank sequence for high-qual clusters
        subProcTup = set([executor.submit(outputRefFlank, \
                                          highQualArray, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 6. Annotate clusters
        subProcTup = set([executor.submit(annotateCluster, \
                                          highQualArray, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
    
    return allCltData