import os
import numpy as np
import subprocess
from collections import OrderedDict, Counter

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

        # 1. Primary assembling
        prefix = "tmp_assm/{}_{}".format(cltView[i].tid, cltView[i].idx)
        cmd = "wtdbg2 -p 5 -k 15 -l 256 -e {0} -S 1 -A --rescue-low-cov-edges --node-len {1} --ctg-min-length {1} " \
            "--ctg-min-nodes 1 -q -t {2} -i {3}.fa -fo {3}".format(minEdge, nodeLen, numThread, prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
        if os.path.isfile(f"{prefix}.ctg.lay.gz") == False:
            continue
        
        cmd = "wtpoa-cns -q -c 1 -t {0} -i {1}.ctg.lay.gz -fo {1}_assm.fa".format(numThread, prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
        if os.path.getsize(f"{prefix}_assm.fa") == 0:
            continue

        # 2. Polishing
        cmd = "minimap2 -aY {0}_assm.fa {0}.fa | samtools sort | samtools view -bhS -o {0}_RawToAssm.bam".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
        if os.path.isfile("f{prefix}_RawToAssm.bam") == False:
            os.rename(f"{prefix}_assm.fa", "f{prefix}_assembled.fa")
            continue

        cmd = "samtools consensus --ff 3332 --homopoly-score 0.1 --low-MQ 10 --scale-MQ 1 --het-scale 0 --P-indel 2e-4 " \
              "{0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 3. Remove N-bases
        cmd = (
            "sed 's/N//g' {0}_assembled.fa | "
            "awk '{{if($1~/^>/){{ctg=$1; a[ctg]=\"\"}} else{{a[ctg]=a[ctg]\"\"$1}}}} "
            "END{{for(ctg in a){{if(a[ctg]!=\"\"){{print ctg; print a[ctg]}}}}}}' > {0}_tmp.fa && "
            "mv {0}_tmp.fa {0}_assembled.fa"
        ).format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 4. Recalibration
        recalibration(prefix, cmdArgs)


cdef int checkLeftSide(int[::1] queryArr, int[::1] refArr,  int i, int leftMost):
    while i >= 0:
        queryPos = queryArr[i]
        refPos = refArr[i]
        if (refPos > 0) and (refPos < leftMost):
            return -1
        if (queryPos < 0) or (refPos < 0):
            i -= 1
            continue
        return queryPos

    return -1


cdef int checkRightSide(int[::1] queryArr, int[::1] refArr, int numPairs,  int i, int rightMost):
    while i < numPairs:
        queryPos = queryArr[i]
        refPos = refArr[i]
        if (refPos > 0) and (refPos > rightMost):
            return -1
        if (queryPos < 0) or (refPos < 0):
            i += 1
            continue
        return queryPos

    return -1


cdef getQueryPolymerLens(BamFile inputBam, int tid, int[::1] queryArr, int[::1] refArr, object polymerRegions):
    cdef Iterator iterator = Iterator(inputBam, tid)
    cdef int queryPos, refPos
    cdef int refStart, refEnd
    cdef int queryStart, queryEnd
    cdef int i, retValue, readLen, numPairs
    
    while True:
        # 1. Load read
        retValue = iterator.cnext1()
        if retValue < 0:
            del iterator
            return

        # 2. Get aligned-pairs
        numPairs = getAlignedPairs(iterator.bamRcord, &queryArr[0], &refArr[0])
        readLen = iterator.bamRcord.core.l_qseq
        refEnd = -2

        # 3. Collect homo-polymer lengths
        for i in range(numPairs):
            queryPos = queryArr[i]
            refPos = refArr[i]

            # 4. Read cover homo-polymer start
            if refPos in polymerRegions:
                queryStart = queryPos
                refStart = refPos
                refEnd = polymerRegions[refPos][0]
                if queryStart < 0:
                    queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, refEnd)
                    if queryPos < 0:
                        refEnd = -2
                    else:
                        queryStart = queryPos
                    continue
                else:
                    queryPos = checkLeftSide(queryArr, refArr, i-1, 0)
                    if queryPos >= 0:
                        queryStart = queryPos + 1            
            
            # 5. Read cover homo-polymer end
            elif refPos == refEnd:
                queryEnd = queryPos
                if queryEnd < 0:
                    queryPos = checkLeftSide(queryArr, refArr, i-1, refStart)
                    if queryPos < 0:
                        refEnd = -2
                    else:
                        queryEnd = queryPos
                    continue
                else:
                    queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, readLen)
                    if queryPos >= 0:
                        queryEnd = queryPos - 1
                
                # 6. Collect query homo-polymer length
                polymerRegions[refStart].append(queryEnd - queryStart + 1)
                refEnd = -2


cdef object getPolymerRegions(str refSeq, int refLen, int minPolymerLen):
    cdef int start = 0
    cdef int end = 0
    cdef object regions = OrderedDict()

    # 1. Find all homo-polymer region
    while end < refLen:
        if refSeq[end] == refSeq[start]:
            end += 1
            continue
        if end - start >= minPolymerLen:
            regions[start] = [end-1]
        
        start = end
        end += 1
    
    # 2. Check final homo-polymer
    if (refSeq[end-1] == refSeq[start]) and (end - start >= minPolymerLen):
        regions[start] = [end-1]
        
    return regions


cdef recalibration(str prefix, object cmdArgs):
    if cmdArgs.recalibration == False:
        return

    # 1. Map raw reads to polished seqs
    os.rename(f"{prefix}_assembled.fa", f"{prefix}_polished.fa")
    cmd = "minimap2 -aY {0}_polished.fa {0}.fa | " \
          "samtools sort | samtools view -bhS -F 3332 -o {0}_RawToPolish.bam && " \
          "samtools index {0}_RawToPolish.bam".format(prefix)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
    if not os.path.exists(f"{prefix}_RawToPolish.bam.bai"):
        os.rename(f"{prefix}_polished.fa", f"{prefix}_assembled.fa")
        return
    
    # 2. Load files
    cdef bytes inputFn = f"{prefix}_polished.fa".encode("utf-8")
    cdef faidx_t *inputFa = fai_load(inputFn)
    cdef BamFile inputBam = BamFile(f"{prefix}_RawToPolish.bam", "rb")
    cdef object outputFa = open(f"{prefix}_assembled.fa", "w")

    # 3. Recalibrate
    cdef object queryArr = np.zeros(1, dtype=np.int32)
    cdef object refArr = np.zeros(1, dtype=np.int32)
    cdef int tid, refLen, newCapacity, numRef = faidx_nseq(inputFa)
    cdef int prevEnd, refStart, refEnd, poymerLen
    cdef const char *cRefName
    cdef char *cRefSeq
    for tid in range(numRef):
        # 4. Load sequence
        cRefName = faidx_iseq(inputFa, tid)
        refLen = faidx_seq_len(inputFa, cRefName)
        cRefSeq = faidx_fetch_seq(inputFa, cRefName, 0, refLen, &refLen)
        refSeq = cRefSeq.decode("utf-8")
        
        # 5. Define homopolymer regions
        polymerRegions = getPolymerRegions(refSeq, refLen, cmdArgs.minPolymerLen)

        # 6. Expand array capacity
        newCapacity = int(2 * refLen)
        if queryArr.shape[0] < newCapacity:
            queryArr.resize((newCapacity,), refcheck=False)
            refArr.resize((newCapacity,), refcheck=False)

        # 7. Output original sequence if no homopolymer
        refName = cRefName.decode('utf-8')
        if len(polymerRegions) == 0:
            outputFa.write(">" + refName + "\n" + refSeq + "\n")
            continue
        
        # 8. Collect query homo-polymer length
        getQueryPolymerLens(inputBam, tid, queryArr, refArr, polymerRegions)

        # 9. Output recalibrated sequence
        prevEnd = 0
        outputFa.write(">" + refName + "\n")
        for refStart in polymerRegions:
            refEnd = polymerRegions[refStart][0]
            outputFa.write(refSeq[prevEnd:refStart])
            
            if len(polymerRegions[refStart]) > 1:
                poymerLen = Counter(polymerRegions[refStart][1:]).most_common(1)[0][0]
            else:
                poymerLen = refEnd - refStart + 1
            
            outputFa.write(refSeq[refStart] * poymerLen)
            prevEnd = refEnd + 1
        
        if prevEnd < len(refSeq):
            outputFa.write(refSeq[prevEnd:])
        outputFa.write("\n")

    # 10. Close files
    fai_destroy(inputFa)
    outputFa.close()
    inputBam.close()