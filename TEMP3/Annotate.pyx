import os
import numpy as np
from subprocess import Popen, DEVNULL

##################
### Data Types ###
##################
AnnoDt = np.dtype([
    ('idx',         np.int32),
    ('refStart',    np.int32),
    ('refEnd',      np.int32),
    ('tid',         np.int32),
    ('strand',      np.uint8),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
])


#########################
### Annotate Assembly ###
#########################
cdef annotateAssm(Cluster[::1] cltArray, int startIdx, int endIdx, object cmdArgs):
    cdef int i

    for i in range(startIdx, endIdx):
        mapFlankToAssm(cltArray[i].tid, cltArray[i].idx)
        outputInsSeq(&cltArray[i])
        mapInsToTE(cltArray[i].tid, cltArray[i].idx, cmdArgs)
        

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


cdef mapInsToTE(int tid, int idx, object cmdArgs):
    cdef int exitCode
    cdef str queryFn = "tmp_anno/{}_{}_insertion.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_InsToTE.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(cmdArgs.teFn, queryFn, outFn)
    if os.path.isfile(queryFn) == False:
        return

    process = Popen([cmd], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitCode = process.wait()
    if exitCode != 0:
        raise Exception("Error: minimap2 failed for {}".format(queryFn))


##########################
### Annotate Insertion ###
##########################
# cdef annotateIns(Cluster[::1] cltArray, int startIdx, int endIdx):
#     cdef int i, returnValue, numAnno=0, maxNum=9900
#     cdef object annoArray = np.zeros(10000, dtype=AnnoDt)
#     cdef Anno[::1] annoArrayView = annoArray

#     for i in range(startIdx, endIdx):
#         if numAnno >= maxNum:
#             maxNum = annoArrayView.shape[0] + 10000
#             annoArray.resize((maxNum,), refcheck=False)
#             annoArrayView = annoArray
#             maxNum -= 100

#         returnValue = fillAnnoArray(&cltArray[i], &annoArrayView[0])
#         numAnno += returnValue
    
#     annoArray.resize((numAnno,), refcheck=False)
#     annoArray.sort(order=['idx', 'queryStart', 'queryEnd'])


###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltArray.shape[0]:
        endIdx = cltArray.shape[0]
    
    annotateAssm(cltArray, startIdx, endIdx, cmdArgs)