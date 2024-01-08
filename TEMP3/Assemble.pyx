import os
import numpy as np
from subprocess import Popen, DEVNULL

######################
### Local Assembly ###
######################
cdef object getAssembleArray(dict chromCltData):
    cdef int i, tid, numClt = 0, maxNumClt = 1900
    cdef object assembleArray = np.zeros((2000, 2), dtype=np.int32)
    cdef int[:, :] arrayView = assembleArray
    cdef Cluster[::1] cltArray

    for tid in chromCltData.keys():
        cltArray = chromCltData[tid][0]
        if cltArray.shape[0] == 0:
            continue
        
        if numClt >= maxNumClt:
            maxNumClt = assembleArray.shape[0] + 2000
            assembleArray.resize((maxNumClt, 2), refcheck=False)
            arrayView = assembleArray
            maxNumClt -= 100
        
        for i in range(cltArray.shape[0]):
            if isLowQualClt(&cltArray[i]) or isSomaClt(&cltArray[i]):
                continue
            
            arrayView[numClt, 0] = tid
            arrayView[numClt, 1] = i
            numClt += 1
    
    assembleArray.resize((numClt, 2), refcheck=False)
    return assembleArray


cpdef wtdbg2Assemble(int[:, :] assembleArray, int start, int taskSize, int numThread):
    cdef int i, exitCode, end
    cdef str cmd, prefix
    cdef object subProcess

    end = start + taskSize
    if end > assembleArray.shape[0]:
        end = assembleArray.shape[0]

    for i in range(start, end):
        prefix = "tmp_assm/tmp.{}_{}".format(assembleArray[i][0], assembleArray[i][1])

        cmd = "wtdbg2 -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 " \
              "--ctg-min-nodes 1 -q -t {} -i {}.fa -fo {}".format(numThread, prefix, prefix)
        subProcess = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitCode = subProcess.wait()

        if os.path.isfile("{}.ctg.lay.gz".format(prefix)) == False:
            continue
        
        cmd = "wtpoa-cns -q -t {} -i {}.ctg.lay.gz -fo {}_assembled.fa".format(numThread, prefix, prefix)
        subProcess = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitCode = subProcess.wait()