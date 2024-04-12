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
    ('flag',        np.uint32),
])


#########################
### Annotate Assembly ###
#########################
cdef annotateAssm(Cluster[::1] cltView, int startIdx, int endIdx, object cmdArgs):
    cdef int i
    for i in range(startIdx, endIdx):
        mapFlankToAssm(cltView[i].tid, cltView[i].idx)
        extractIns(&cltView[i])
        mapInsToTE(cltView[i].tid, cltView[i].idx, cmdArgs)
        

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


##########################
### Annotate Insertion ###
##########################
cdef object getClassArray(object cmdArgs):
    cdef int numTE=0, maxNum=190
    cdef object classArray = np.zeros(200, dtype=np.uint32)
    cdef uint32_t[::1] classView = classArray

    for l in open(cmdArgs.classFn, "r"):
        l = l.strip().split()
        if len(l) < 2:
            continue
        
        if numTE >= maxNum:
            maxNum = classView.shape[0] + 200
            classArray.resize((maxNum,), refcheck=False)
            classView = classArray
            maxNum -= 10
        
        if l[1] == "DNA":
            classView[numTE] = CLT_DNA
        elif l[1] == "LTR":
            classView[numTE] = CLT_LTR
        elif l[1] == "LINE":
            classView[numTE] = CLT_LINE
        elif l[1] == "SINE":
            classView[numTE] = CLT_SINE
        elif l[1] == "RETROPOSON":
            classView[numTE] = CLT_RETROPOSON
        
        numTE += 1
    
    classArray.resize((numTE,), refcheck=False)
    return classArray


cdef object annotateIns(Cluster[::1] cltView, int startIdx, int endIdx, object cmdArgs):
    cdef int i, numTmp, numAnno=0, maxNum=2900
    cdef object annoArray = np.zeros(3000, dtype=AnnoDt)
    cdef object classArray = getClassArray(cmdArgs)
    cdef Anno[::1] annoView = annoArray
    cdef uint32_t[::1] classView = classArray
    cdef str bamFn, assmFn

    for i in range(startIdx, endIdx):
        assmFn = "tmp_assm/{}_{}.ctg.lay.gz".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(assmFn) != 0:
                cltView[i].flag |= CLT_ASSEMBLED
        
        bamFn = "tmp_anno/{}_{}_InsToTE.bam".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(bamFn) == False:
            continue

        if numAnno >= maxNum:
            maxNum = annoView.shape[0] + 3000
            annoArray.resize((maxNum,), refcheck=False)
            annoView = annoArray
            maxNum -= 100

        numTmp = fillAnnoArray(&cltView[i], &annoView[numAnno], i)
        mapTsdToLocal(cltView[i].tid, cltView[i].idx)
        bamFn = "tmp_anno/{}_{}_TsdToLocal.bam".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(bamFn) == True:
            annoTsd(&cltView[i], &annoView[numAnno], numTmp)

        setInsStruc(&cltView[i], &annoView[numAnno], numTmp, &classView[0])
        numAnno += numTmp
    
    annoArray.resize((numAnno,), refcheck=False)
    annoArray.sort(order=['idx', 'queryStart', 'queryEnd'])
    return annoArray


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


###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotateCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]
    
    annotateAssm(cltView, startIdx, endIdx, cmdArgs)
    cdef object annoArray = annotateIns(cltView, startIdx, endIdx, cmdArgs)
    
    # Output formated cluster and annotation records
    cdef Anno[::1] annoView = annoArray
    cdef bytes teFn = cmdArgs.teFn.encode()
    cdef bytes refFn = cmdArgs.refFn.encode()
    outputAnno(&annoView[0], annoView.shape[0], startIdx, teFn)
    outputClt(&cltView[0], startIdx, endIdx, refFn, teFn)

    # Output cltArray and annoArray for post-anno-filtering
    cdef object cltArray = np.asarray(cltView)
    annoArray.tofile('tmp_anno/{}_annoArray.dat'.format(startIdx))
    cltArray[startIdx:endIdx].tofile('tmp_anno/{}_cltArray.dat'.format(startIdx))