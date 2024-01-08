import numpy as np
from .Cluster import buildCluster
from .Assemble import wtdbg2Assemble
from .FileIO import outputSomaCltSeqs
from concurrent.futures import ProcessPoolExecutor, as_completed

cdef dict getBackgroundInfo(str genomeBamFilePath, int numThread):
    
    cdef BamFile genomeBamFile = BamFile(genomeBamFilePath, "rb", numThread)
    cdef int i, tid=0, maxChromLen=0

    for i in range(genomeBamFile.header.n_targets):
        if maxChromLen < sam_hdr_tid2len(genomeBamFile.header, i):
            maxChromLen = sam_hdr_tid2len(genomeBamFile.header, i)
            tid = i
    
    cdef Iterator iterator = Iterator(genomeBamFile, tid)
    cdef int numAln = 0, maxNumAln = 499980
    cdef object templatei = np.zeros(500000, dtype=np.int32)
    cdef object templatef = np.zeros(500000, dtype=np.float32)
    cdef object readLenArray = np.zeros(500000, dtype=np.int32)
    cdef object divArray = np.zeros(500000, dtype=np.float32)
    cdef float[::1] divArrayView = divArray
    cdef int[::1] readLenArrayView = readLenArray
    cdef int returnValue, alnLen
    cdef int64_t sumAlnLen = 0
    cdef float divergence
    cdef dict bgInfo = {}

    while True:
        returnValue = iterator.cnext1()
        if returnValue < 0:
            bgInfo["numChrom"] = genomeBamFile.header.n_targets
            bgInfo["bgDiv"] = np.mean(divArray[:numAln])
            bgInfo["bgDepth"] = float(sumAlnLen) / maxChromLen
            bgInfo["bgReadLen"] = np.median(readLenArray[:numAln])

            genomeBamFile.close(); del iterator; del genomeBamFile; return bgInfo

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numAln > maxNumAln:
            divArray = np.concatenate((divArray, templatef))
            readLenArray = np.concatenate((readLenArray, templatei))
            divArrayView = divArray
            readLenArrayView = readLenArray
            maxNumAln = divArray.shape[0] - 20
            
        getMapLenAndDiv(&alnLen, &divergence, iterator.bamRcord)
        sumAlnLen += alnLen
        divArrayView[numAln] = divergence
        readLenArrayView[numAln] = iterator.bamRcord.core.l_qseq
        numAln += 1


cdef tuple divideTasks(int numTasks, int poolSize):
    cdef int i, taskSize = int(numTasks / poolSize)
    cdef list startList = []

    if taskSize == 0:
        taskSize = 1
    
    for i in range(0, numTasks, taskSize):
        startList.append(i)
    
    return taskSize, startList


cpdef dict runInParallel(object cmdArgs):
    
    cdef set subProcTup
    cdef dict allCltData = {}, chromCltData, bgInfo
    cdef object subProc, assembleArray
    cdef int start, taskSize
    cdef list startList

    # 1. Get Background Info
    bgInfo = getBackgroundInfo(cmdArgs.genomeBamFilePath, cmdArgs.numThread)
    print("bg Divergence: {}\nbg coverage: {}\nbg readlen: {}" \
          "".format(bgInfo["bgDiv"], bgInfo["bgDepth"], bgInfo["bgReadLen"]))

    with ProcessPoolExecutor(max_workers=cmdArgs.numProcess) as executor:
        # 2. Get Background Info
        subProcTup = set([executor.submit(buildCluster, \
                                          bgInfo["bgDiv"], \
                                          bgInfo["bgDepth"], \
                                          bgInfo["bgReadLen"], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
                                       
        for subProc in as_completed(subProcTup):
            chromCltData = subProc.result()
            allCltData = {**allCltData, **chromCltData}
        
        # 3. Local Assembly
        assembleArray = getAssembleArray(allCltData)
        taskSize, startList = divideTasks(assembleArray.shape[0], cmdArgs.numProcess)
        subProcTup = set([executor.submit(wtdbg2Assemble, \
                                          assembleArray, \
                                          start, \
                                          taskSize, \
                                          cmdArgs.numThread) for start in startList])

        for subProc in as_completed(subProcTup):
            returnValue = subProc.result()
        
        # 4. Output sequence for clusters without assembly
        subProcTup = set([executor.submit(outputSomaCltSeqs, \
                                          allCltData[tid][0], \
                                          allCltData[tid][1], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
                                       
        for subProc in as_completed(subProcTup):
            returnValue = subProc.result()
    
    return allCltData