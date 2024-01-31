import os
import numpy as np
from subprocess import Popen, DEVNULL

##################
### Data Types ###
##################
AnnoDt = np.dtype([
    ('idx',         np.int32),
    ('cltTid',      np.int32),
    ('cltIdx',      np.int32),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
    ('strand',      np.uint8),
    ('tid',         np.int32),
    ('refStart',    np.int32),
    ('refEnd',      np.int32),
])


#########################
### Annotate Assembly ###
#########################
cdef annotateAssm(Cluster[::1] cltArray, int startIdx, int endIdx, object cmdArgs):
    cdef int i

    for i in range(startIdx, endIdx):
        mapFlankToAssm(cltArray[i].tid, cltArray[i].idx)
        extractIns(&cltArray[i])
        mapInsToTE(cltArray[i].tid, cltArray[i].idx, cmdArgs)
        mapTsdToLocal(cltArray[i].tid, cltArray[i].idx)
        

cdef mapFlankToAssm(int tid, int idx):
    cdef int exitCode
    cdef str targetFn = "tmp_assm/{}_{}_assembled.fa".format(tid, idx)
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


cdef mapTsdToLocal(int tid, int idx):
    cdef int exitCode
    cdef str targetFn = "tmp_anno/{}_{}_local.fa".format(tid, idx)
    cdef str queryFn = "tmp_anno/{}_{}_tsd.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_TsdToLocal.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(targetFn, queryFn, outFn)
    if os.path.isfile(queryFn) == False:
        return

    process = Popen([cmd], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitCode = process.wait()
    if exitCode != 0:
        raise Exception("Error: minimap2 failed for {}".format(queryFn))


##########################
### Annotate Insertion ###
##########################
cdef object annotateIns(Cluster[::1] cltArray, int startIdx, int endIdx):
    cdef int i, retValue, numAnno=0, maxNum=9900
    cdef object annoArray = np.zeros(10000, dtype=AnnoDt)
    cdef Anno[::1] annoArrayView = annoArray
    cdef str bamFn, assmFn

    for i in range(startIdx, endIdx):
        assmFn = "tmp_assm/{}_{}.ctg.lay.gz".format(cltArray[i].tid, cltArray[i].idx)
        if os.path.isfile(assmFn) != 0:
                cltArray[i].flag |= CLT_ASSEMBLED
        
        bamFn = "tmp_anno/{}_{}_InsToTE.bam".format(cltArray[i].tid, cltArray[i].idx)
        if os.path.isfile(bamFn) == False:
            continue

        if numAnno >= maxNum:
            maxNum = annoArrayView.shape[0] + 10000
            annoArray.resize((maxNum,), refcheck=False)
            annoArrayView = annoArray
            maxNum -= 100

        retValue = fillAnnoArray(&cltArray[i], &annoArrayView[numAnno], i)
        annoTsd(&cltArray[i])
        numAnno += retValue
    
    annoArray.resize((numAnno,), refcheck=False)
    annoArray.sort(order=['idx', 'queryStart', 'queryEnd'])
    annoArrayView = annoArray

    cdef bytes annoFn = 'tmp_anno/{}_anno.txt'.format(startIdx).encode()
    outPutAnno(&annoArrayView[0], numAnno, annoFn)

    ##
    return annoArray
    ##


###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltArray.shape[0]:
        endIdx = cltArray.shape[0]
    
    annotateAssm(cltArray, startIdx, endIdx, cmdArgs)

    ##
    annoArray = annotateIns(cltArray, startIdx, endIdx)
    # np.savetxt('tmp_anno_{}.txt'.format(startIdx), annoArray, fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d')
    np.savetxt('tmp_clt_{}.txt'.format(startIdx), cltArray[startIdx:endIdx-1])

    rev = 0; black = 0; assm = 0; left = 0; right = 0; diff = 0; same = 0; te = 0; polya = 0; tsd = 0
    for i in range(startIdx, endIdx):
        if (cltArray[i].flag & 1) != 0:
            rev += 1
        if (cltArray[i].flag & 2) != 0:
            black += 1
        if (cltArray[i].flag & 4) != 0:
            assm += 1
        if (cltArray[i].flag & 8) != 0:
            left += 1
        if (cltArray[i].flag & 16) != 0:
            right += 1
        if (cltArray[i].flag & 32) != 0:
            diff += 1
        if (cltArray[i].flag & 64) != 0:
            same += 1
        if (cltArray[i].flag & 128) != 0:
            te += 1
        if (cltArray[i].flag & 256) != 0:
            polya += 1
        if (cltArray[i].flag & 512) != 0:
            tsd += 1
    
    print("rev:{}\tblack:{}\tassm:{}\tleft:{}\tright:{}\tdiff:{}\tsame:{}\tte:{}\tpolya:{}\ttsd:{}".format(rev, black, assm, left, right, diff, same, te, polya, tsd))
    ##