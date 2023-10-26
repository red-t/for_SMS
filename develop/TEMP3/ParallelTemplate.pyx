import numpy as np
from .Cluster import buildCluster
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
    cdef object deArray = np.zeros(500000, dtype=np.float32)
    cdef float[::1] divArrayView = divArray
    cdef float[::1] deArrayView = deArray
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
            bgInfo["bgDe"] = np.mean(deArray[:numAln])
            bgInfo["bgDepth"] = float(sumAlnLen) / maxChromLen
            bgInfo["bgReadLen"] = np.median(readLenArray[:numAln])

            genomeBamFile.close(); del iterator; del genomeBamFile; return bgInfo

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numAln > maxNumAln:
            divArray = np.concatenate((divArray, templatef))
            deArray = np.concatenate((deArray, templatef))
            readLenArray = np.concatenate((readLenArray, templatei))
            divArrayView = divArray
            deArrayView = deArray
            readLenArrayView = readLenArray
            maxNumAln = deArray.shape[0] - 20
            
        getMapLenAndDiv(&alnLen, &divergence, iterator.bamRcord)
        sumAlnLen += alnLen
        divArrayView[numAln] = divergence
        deArrayView[numAln] = getDe(iterator.bamRcord)
        readLenArrayView[numAln] = iterator.bamRcord.core.l_qseq
        numAln += 1


cpdef dict buildClusterParallel(object args):
    '''call buildCluster in multi-process way'''
    cdef set futures
    cdef dict result = {}, returnValue
    cdef dict bgInfo = getBackgroundInfo(args.genomeBamFilePath, args.numThread)

    print("bg Divergence: {}\nbg de: {}\nbg coverage: {}\nbg readlen: {}".format(bgInfo["bgDiv"], bgInfo["bgDe"], bgInfo["bgDepth"], bgInfo["bgReadLen"]))

    with ProcessPoolExecutor(max_workers=args.numProcess) as executor:
        futures = set([executor.submit(buildCluster,
                                       args.genomeBamFilePath,
                                       args.repeatPath,
                                       args.gapPath,
                                       args.referenceTe,
                                       args.numThread,
                                       tid,
                                       args.minSegLen,
                                       args.maxDistance,
                                       bgInfo["bgDiv"],
                                       bgInfo["bgDe"],
                                       bgInfo["bgDepth"],
                                       bgInfo["bgReadLen"]) for tid in range(bgInfo["numChrom"])])
                                       
        for future in as_completed(futures):
            returnValue = future.result()
            result = {**result, **returnValue}
    
    return result