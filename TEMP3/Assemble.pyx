import os
import numpy as np
import subprocess

######################
### Local Assembly ###
######################
cdef int getMinEdge(int numSegRaw):
    if numSegRaw < 5:
        return 1
    elif numSegRaw < 10:
        return 2
    elif numSegRaw < 40:
        return 3
    elif numSegRaw < 70:
        return 4
    else:
        return 5


cpdef assembleCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs):
    cdef int numThread = cmdArgs.numThread
    cdef int nodeLen = cmdArgs.nodeLen
    cdef int i, endIdx, minEdge
    cdef str cmd, prefix

    endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]

    for i in range(startIdx, endIdx):
        if isSomaClt(&cltView[i]):
            continue

        minEdge = getMinEdge(cltView[i].numSegRaw)
        if cmdArgs.minEdge > 0:
            minEdge = cmdArgs.minEdge

        # Primary assembling    
        prefix = "tmp_assm/{}_{}".format(cltView[i].tid, cltView[i].idx)
        cmd = "wtdbg2 -p 5 -k 15 -l 256 -e {0} -S 1 -A --rescue-low-cov-edges --node-len {1} --ctg-min-length {1} " \
            "--ctg-min-nodes 1 -q -t {2} -i {3}.fa -fo {3}".format(minEdge, nodeLen, numThread, prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        if os.path.isfile("{}.ctg.lay.gz".format(prefix)) == False:
            continue
        
        cmd = "wtpoa-cns -q -c 1 -t {0} -i {1}.ctg.lay.gz -fo {1}_assm.fa".format(numThread, prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # Polishing
        cmd = "minimap2 -aY {0}_assm.fa {0}.fa | samtools sort | samtools view -bhS -o {0}_RawToAssm.bam".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        if os.path.isfile("{}_RawToAssm.bam".format(prefix)) == False:
            os.rename("{}_assm.fa".format(prefix), "{}_assembled.fa".format(prefix))
            continue

        cmd = "samtools consensus --homopoly-score 0.1 --low-MQ 10 --scale-MQ 1 --het-scale 0 --P-indel 2e-4 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # Remove N-bases
        cmd = (
            "sed 's/N//g' {0}_assembled.fa | "
            "awk '{{if($1~/^>/){{ctg=$1; a[ctg]=\"\"}} else{{a[ctg]=a[ctg]\"\"$1}}}} "
            "END{{for(ctg in a){{if(a[ctg]!=\"\"){{print ctg; print a[ctg]}}}}}}' > {0}_tmp.fa && "
            "mv {0}_tmp.fa {0}_assembled.fa"
        ).format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')