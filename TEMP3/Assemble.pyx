import os
import numpy as np
import subprocess
from collections import OrderedDict, Counter, defaultdict
from cpython cimport PyBytes_FromStringAndSize

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

        # 2. First round polishing
        cmd = "minimap2 -aY {0}_assm.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
              "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 3. Second round polishing
        cmd = "minimap2 -aY {0}_assembled.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
              "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 4. Recalibration
        recalibration(prefix, cmdArgs)

        # 5. Remove N-bases
        cmd = (
            "sed 's/N//g' {0}_assembled.fa | "
            "awk '{{if($1~/^>/){{ctg=$1; a[ctg]=\"\"}} else{{a[ctg]=a[ctg]\"\"$1}}}} "
            "END{{for(ctg in a){{if(a[ctg]!=\"\"){{print ctg; print a[ctg]}}}}}}' > {0}_tmp.fa && "
            "mv {0}_tmp.fa {0}_assembled.fa"
        ).format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


#####################
### Recalibration ###
#####################
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
    cdef int refStart, refEnd
    cdef const char *cRefName
    cdef char *cRefSeq

    for tid in range(numRef):
        # 4. Load sequence
        cRefName = faidx_iseq(inputFa, tid)
        refLen = faidx_seq_len(inputFa, cRefName)
        cRefSeq = faidx_fetch_seq(inputFa, cRefName, 0, refLen, &refLen)
        refSeq = cRefSeq.decode('utf-8')
        
        # 5. Define homopolymer regions
        polymerRegions = getPolymerRegions(refSeq, refLen)

        # 6. Expand array capacity
        newCapacity = int(2 * refLen)
        if queryArr.shape[0] < newCapacity:
            queryArr.resize((newCapacity,), refcheck=False)
            refArr.resize((newCapacity,), refcheck=False)
        
        # 8. Collect query homo-polymers
        iterator = Iterator(inputBam, tid)
        queryPolymers = getQueryPolymers(iterator, queryArr, refArr, polymerRegions)
        
        # 9. Output recalibrated sequence
        outputFa.write(">" + cRefName.decode('utf-8') + "\n")
        for refStart in polymerRegions:
            # Ignore long homopolymer, un-covered region
            refEnd = polymerRegions[refStart]
            if (refEnd - refStart + 1 >= 20) or (refStart not in queryPolymers):
                outputFa.write(refSeq[refStart:refEnd+1])
                continue

            outputFa.write(getRecalibratedSeq(queryPolymers, refSeq, refStart, refEnd))
        
        outputFa.write("\n")
        del iterator

    # 10. Close files
    fai_destroy(inputFa)
    outputFa.close()
    inputBam.close()


cdef object getPolymerRegions(str refSeq, int refLen):
    cdef int start = 0
    cdef int end = 1
    cdef object regions = OrderedDict()

    while end < refLen:
        if refSeq[end] == refSeq[start]:
            end += 1
            continue
        
        regions[start] = end - 1
        start = end
        end += 1
    
    # Final homo-polymer
    regions[start] = end - 1
    return regions


cdef object getQueryPolymers(Iterator iterator, int[::1] queryArr, int[::1] refArr, object polymerRegions):
    cdef int refPos, refStart, refEnd
    cdef int queryPos, queryStart, queryEnd
    cdef int i, retValue, readLen, numPairs
    cdef object queryPolymers = OrderedDict()
    cdef str readSeq, querySeq
    
    while True:
        # 1. Load read
        retValue = iterator.cnext1()
        if retValue < 0:
            return queryPolymers

        # 2. Get aligned-pairs
        readSeq = getReadSeq(iterator)
        readLen = iterator.bamRcord.core.l_qseq
        numPairs = getAlignedPairs(iterator.bamRcord, &queryArr[0], &refArr[0])
        refEnd = -2

        # 3. Collect homo-polymer lengths
        for i in range(numPairs):
            queryPos = queryArr[i]
            refPos = refArr[i]

            # 4. Read cover homo-polymer start
            if refPos in polymerRegions:
                queryStart = queryPos
                refStart = refPos
                refEnd = polymerRegions[refPos]

                # Single base
                if refStart == refEnd:
                    if queryStart < 0:
                        queryStart = -1
                        queryEnd = -2
                        querySeq = ""
                    else:
                        queryStart = checkLeftSide(queryArr, refArr, i-1, 0) + 1
                        queryEnd = queryPos
                        querySeq = readSeq[queryStart:queryEnd+1]
                    
                    if refStart in queryPolymers:
                        queryPolymers[refStart].append((querySeq, queryEnd - queryStart + 1))
                    else:
                        queryPolymers[refStart] = [(querySeq, queryEnd - queryStart + 1)]
                    continue
                
                # Homo-polymer region
                if queryStart < 0:
                    # Correction for boundary deletion
                    queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, refEnd)
                    if queryPos < 0:
                        refEnd = -2
                    else:
                        queryStart = queryPos
                    continue
                else:
                    # Correction for boundary insertion
                    queryPos = checkLeftSide(queryArr, refArr, i-1, 0)
                    if queryPos >= 0:
                        queryStart = queryPos + 1            
            
            # 5. Read cover homo-polymer end
            elif refPos == refEnd:
                queryEnd = queryPos
                if queryEnd < 0:
                    # Correction for boundary deletion
                    queryPos = checkLeftSide(queryArr, refArr, i-1, refStart)
                    if queryPos < 0:
                        refEnd = -2
                    else:
                        queryEnd = queryPos
                    continue
                
                # 6. Collect query homo-polymers
                if refStart in queryPolymers:
                    queryPolymers[refStart].append((readSeq[queryStart:queryEnd+1], queryEnd - queryStart + 1))
                else:
                    queryPolymers[refStart] = [(readSeq[queryStart:queryEnd+1], queryEnd - queryStart + 1)]
                
                refEnd = -2


cdef str getReadSeq(Iterator iterator):
    cdef int i, readLen = iterator.bamRcord.core.l_qseq
    cdef bytes bSeq = PyBytes_FromStringAndSize(NULL, readLen)
    cdef char *cSeq = <char*>bSeq

    for i in range(readLen):
        cSeq[i] = seq_nt16_str[bam_seqi(bam_get_seq(iterator.bamRcord), i)]

    return bSeq.decode('utf-8')


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


cdef int getMostCommonLen(object queryPolymers, int refStart, int refEnd):
    cdef object counter = Counter([p[1] for p in queryPolymers[refStart]])
    cdef int refLen = refEnd - refStart + 1
    cdef int refCount = counter[refLen]
    cdef int altLen, altCount

    for altLen, altCount in counter.most_common(2):
        if altLen != refLen:
            break
    
    if altCount < 0.8 * refCount:
        return refLen
    else:
        return altLen


cdef str getRecalibratedSeq(object queryPolymers, str refSeq, int refStart, int refEnd):
    # 1. Get most common length
    cdef int mostCommonLen = getMostCommonLen(queryPolymers, refStart, refEnd)
    if mostCommonLen == 0:
        return ''

    # 2. Initialize counter for each position
    cdef list positionCounters = [defaultdict(int) for _ in range(mostCommonLen)]
    cdef object polymer
    cdef int i, queryLen
    cdef str querySeq, base
    
    # 3. Count for each position
    for polymer in queryPolymers[refStart]:
        queryLen = polymer[1]
        if queryLen != mostCommonLen:
            continue
        
        querySeq = polymer[0]
        for i in range(mostCommonLen):
            base = querySeq[i]
            positionCounters[i][base] += 1
    
    # 4. Get most common base for each position
    cdef list result = []
    cdef object counter
    if len(positionCounters[0]) == 0:
        return refSeq[refStart:refEnd+1]
    else:
        for counter in positionCounters:
            base = max(counter, key=counter.get)
            result.append(base)

        return ''.join(result)