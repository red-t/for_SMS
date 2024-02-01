import os
import numpy as np
from subprocess import Popen, DEVNULL

######################
### Local Assembly ###
######################
cpdef assembleCluster(Cluster[::1] cltView, int startIdx, int taskSize, int numThread):
    cdef int i, exitCode, endIdx
    cdef str cmd, prefix
    cdef object subProcess

    endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]

    for i in range(startIdx, endIdx):
        if isSomaClt(&cltView[i]):
            continue

        prefix = "tmp_assm/{}_{}".format(cltView[i].tid, cltView[i].idx)
        cmd = "wtdbg2 -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 " \
              "--ctg-min-nodes 1 -q -t {} -i {}.fa -fo {}".format(numThread, prefix, prefix)
        subProcess = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitCode = subProcess.wait()

        if os.path.isfile("{}.ctg.lay.gz".format(prefix)) == False:
            continue
        
        cmd = "wtpoa-cns -q -t {} -i {}.ctg.lay.gz -fo {}_assembled.fa".format(numThread, prefix, prefix)
        subProcess = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitCode = subProcess.wait()