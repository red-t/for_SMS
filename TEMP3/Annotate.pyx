from subprocess import Popen, DEVNULL

#########################
### Annotate Assembly ###
#########################
cdef annotateAssm(Cluster[::1] cltArray, int startIdx, int endIdx):
    cdef int i

    for i in range(startIdx, endIdx):
        mapFlankToAssm(cltArray[i].tid, cltArray[i].idx)
        outputInsSeq(&cltArray[i])
        

cdef mapFlankToAssm(int tid, int idx):
    cdef int exitCode
    cdef str targetFn = "tmp_assm/tmp.{}_{}_assembled.fa".format(tid, idx)
    cdef str queryFn = "tmp_anno/{}_{}_flank.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_FlankToAssm.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(targetFn, queryFn, outFn)

    process = Popen([cmd], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitCode = process.wait()
    if exitCode != 0:
        raise Exception("Error: minimap2 failed for {}".format(queryFn))


###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltArray.shape[0]:
        endIdx = cltArray.shape[0]
    
    annotateAssm(cltArray, startIdx, endIdx)