import os

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