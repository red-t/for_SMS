import os

########################
### Class Definition ###
########################
cdef int MAX_POS = (1 << 31) - 1
cdef class BamFile:
    def __cinit__(self, str filePath, str mode, int numThread=1, BamFile template=None):
        cdef bytes filePathBytes = os.fsencode(filePath)
        cdef bytes indexPathBytes = os.fsencode(filePath + ".bai")
        cdef bytes modeBytes = mode.encode()

        self.filePath = filePathBytes
        self.indexFilePath = indexPathBytes
        self.mode = modeBytes
        self.numThread = numThread

        self.openFile(template)
    
    cdef void openFile(self, BamFile template=None):
        self.htsFile = self.openHtsFile()
        if self.htsFile == NULL:
            raise IOError("could not open BAM file `{}`".format(self.filePath.decode()))

        if self.mode == b'rb':
            with nogil:
                self.header = sam_hdr_read(self.htsFile)
            if self.header == NULL:
                raise IOError("file does not have a valid header, is it BAM format?")
            
            if os.path.exists(self.indexFilePath.decode()):
                with nogil:
                    self.index = sam_index_load2(self.htsFile, self.filePath, self.indexFilePath)
                if not self.index:
                    raise IOError('unable to open index file `{}`'.format(self.indexFilePath.decode()))

        elif (self.mode == b'wb') or (self.mode == b'wF'):
            if template:
                self.header = sam_hdr_dup(template.header)
            else:
                raise ValueError("need a template for copying header")

            if self.mode == b'wb':
                with nogil:
                    sam_hdr_write(self.htsFile, self.header)
    
    cdef htsFile *openHtsFile(self) except? NULL:
        '''open file in 'rb/wb/wF' mode'''
        cdef htsFile *htsFile
        with nogil:
            htsFile = sam_open(self.filePath, self.mode)
            if htsFile != NULL:
                hts_set_threads(htsFile, self.numThread)
            return htsFile
    
    cdef void write(self, bam1_t *bamRecord):
        cdef int returnValue
        with nogil:
            returnValue = sam_write1(self.htsFile, self.header, bamRecord)
        if returnValue < 0:
            raise IOError("sam_write1 failed with error code {}".format(returnValue))
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def close(self):
        if self.htsFile:
            sam_close(self.htsFile)
            self.htsFile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.header:
            sam_hdr_destroy(self.header)
            self.header = NULL

    def __dealloc__(self):
        if self.htsFile:
            sam_close(self.htsFile)
            self.htsFile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.header:
            sam_hdr_destroy(self.header)
            self.header = NULL


cdef class Iterator:
    def __cinit__(self):
        self.bamRcord = bam_init1()
        if self.bamRcord == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self, BamFile bamFile, int tid=-1):
        self.htsFile = bamFile.htsFile
        if tid >= 0:
            with nogil:
                self.iter = sam_itr_queryi(bamFile.index, tid, 0, MAX_POS)
    
    def __dealloc__(self):
        if self.bamRcord:
            bam_destroy1(self.bamRcord)
            self.bamRcord = NULL
        if self.iter:
            sam_itr_destroy(self.iter)
            self.iter = NULL
    
    def __iter__(self):
        return self

    cdef int cnext1(self):
        '''cversion of iterator. returnValue>=0 if success.'''
        cdef int returnValue

        with nogil:
            if self.iter.curr_off == 0:
                if self.iter.n_off > 0:
                    self.offset = self.iter.off[0].u
                else:
                    self.offset = self.iter.curr_off
            else:
                self.offset = self.iter.curr_off

            returnValue = hts_itr_next(self.htsFile.fp.bgzf, self.iter, self.bamRcord, self.htsFile)

        return returnValue
    
    cdef int cnext2(self):
        '''directly read a alignment record'''
        cdef int returnValue

        with nogil:
            returnValue = bam_read1(self.htsFile.fp.bgzf, self.bamRcord)
        
        return returnValue

    cdef int cnext3(self, int64_t offset):
        '''read a alignment record with specified offset'''
        cdef int returnValue

        with nogil:
            bgzf_seek(self.htsFile.fp.bgzf, offset, SEEK_SET)
            returnValue = bam_read1(self.htsFile.fp.bgzf, self.bamRcord)
        
        return returnValue


######################
### Construct Args ###
######################
cdef Args newArgs(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs):
    cdef Args args

    args.tid = tid
    args.bgDiv = bgDiv
    args.bgDepth = bgDepth
    args.bgReadLen = bgReadLen
    args.numThread = cmdArgs.numThread
    args.minSegLen = cmdArgs.minSegLen
    args.maxDistance = cmdArgs.maxDistance
    args.minOverhang = cmdArgs.minOverhang
    return args

#######################
### Construc AiList ###
#######################
cdef AiList* newAiList(str bedFn, const char *chrom):
    cdef bytes bedFnBytes = bedFn.encode()
    cdef AiList *aiList = initAiList()

    readBED(aiList, bedFnBytes, chrom)
    constructAiList(aiList, 20)
    return aiList


###########################
### Segment Sequence IO ###
###########################
cdef ouputAllSegSeqs(Segment[::1] segArray, BamFile genomeBam, Args args):
    cdef str outputFileName = "tmp_build/all_seg_{}.fa".format(args.tid)
    cdef BamFile outputFasta = BamFile(outputFileName, "wF", args.numThread, genomeBam)
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, returnValue

    for i in range(segArray.shape[0]):
        returnValue = iterator.cnext3(segArray[i].fileOffset)
        trimSegment(iterator.bamRcord, destRecord, i, segArray[i].queryStart, segArray[i].queryEnd)
        outputFasta.write(destRecord)
    
    bam_destroy1(destRecord); outputFasta.close(); del outputFasta; del iterator


cdef outputGermCltSeqs(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBam, Args args):
    cdef str outputFileName
    cdef BamFile outputFasta
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, j

    for i in range(cltArray.shape[0]):
        if isLowQualClt(&cltArray[i]) or isSomaClt(&cltArray[i]):
            continue

        outputFileName = "tmp_assm/{}_{}.fa".format(args.tid, i)
        outputFasta = BamFile(outputFileName, "wF", args.numThread, genomeBam)
        for j in range(cltArray[i].startIdx, cltArray[i].endIdx):
            if overhangIsShort(&segArray[j], args.minOverhang):
                continue
            outputSingleSeq(segArray, outputFasta, iterator, destRecord, j)
        
        outputFasta.close()

    bam_destroy1(destRecord); del iterator


cpdef outputSomaCltSeqs(Cluster[::1] cltArray, Segment[::1] segArray, object cmdArgs, int tid):
    cdef int i, j
    cdef Args args
    cdef str outputFileName
    cdef BamFile outputFasta
    cdef BamFile genomeBam = BamFile(cmdArgs.genomeBamFilePath, "rb", cmdArgs.numThread)
    cdef Iterator iterator = Iterator(genomeBam, tid)
    cdef bam1_t *destRecord = bam_init1()

    args.minOverhang = cmdArgs.minOverhang
    for i in range(cltArray.shape[0]):
        if isLowQualClt(&cltArray[i]):
            continue
        
        # Skip successfully assembled clusters
        outputFileName = "tmp_assm/{}_{}_assembled.fa".format(tid, i)
        if os.path.isfile(outputFileName):
            if os.path.getsize(outputFileName) != 0:
                cltArray[i].flag |= CLT_ASSEMBLED
                continue
        
        outputFasta = BamFile(outputFileName, "wF", cmdArgs.numThread, genomeBam)
        j = getOuputSegIdx(&cltArray[i], &segArray[0], args)
        outputSingleSeq(segArray, outputFasta, iterator, destRecord, j)
        outputFasta.close()

    bam_destroy1(destRecord); del iterator
    genomeBam.close(); del genomeBam


cdef outputSingleSeq(Segment[::1]segArray, BamFile outputFasta, Iterator iterator, bam1_t *destRecord, int j, int flankSize=3000):
    cdef int start, end, returnValue

    returnValue = iterator.cnext3(segArray[j].fileOffset)
    setTrimRegion(&segArray[j], &start, &end, flankSize)
    trimSegment(iterator.bamRcord, destRecord, j, start, end)
    outputFasta.write(destRecord)


#########################
### Flank Sequence IO ###
#########################
cpdef outputRefFlank(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltArray.shape[0]:
        endIdx = cltArray.shape[0]
    
    cdef bytes refFn = cmdArgs.refFn.encode('utf-8')
    extractRefFlanks(refFn, &cltArray[0], startIdx, endIdx)