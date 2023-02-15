########################################################
## global variables
# maximum genomic coordinace
# for some reason, using 'int' causes overflow
cdef int MAX_POS = (1 << 31) - 1


cdef class BamFile:
    def __cinit__(self, str filepath, int32_t nthreads=1):
        self.threads = nthreads
        self.filename = filename = encode_filename(filepath)
        self.index_filename = idxfilename = encode_filename(filepath+".bai")
        self._open()
        
    def close(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    cdef _open(self):
        '''open BAM file and index file'''
        # open BAM file for reading
        self.htsfile = self._open_htsfile()
        if self.htsfile == NULL:
            if errno:
                raise IOError(errno, "could not open alignment file `{}`: {}".format(force_str(self.filename), force_str(strerror(errno))))
            else:
                raise ValueError("could not open alignment file `{}`".format(force_str(self.filename)))

        # bam files require a valid header
        # cdef bam_hdr_t  *hdr = NULL
        # with nogil:
        #     hdr = sam_hdr_read(self.htsfile)
        # if hdr == NULL:
        #     raise ValueError("file does not have a valid header, is it BAM format?")
        # self.header = makeAlignmentHeader(hdr)
        
        # open BAM index file to enable random access
        with nogil:
            self.index = sam_index_load2(self.htsfile, self.filename, self.index_filename)
        if not self.index:
            if errno:
                raise IOError(errno, force_str(strerror(errno)))
            else:
                raise IOError('unable to open index file `{}`'.format(force_str(self.index_filename)))
    
    cdef htsFile *_open_htsfile(self) except? NULL:
        '''open BAM file in `rb` mode, return htsFile object if success.'''
        cdef char *cfilename
        cdef char *cmode
        cdef int32_t threads = self.threads - 1
        cdef htsFile *htsfile=NULL

        cmode = mode = force_bytes('rb')
        if isinstance(self.filename, bytes):
            cfilename = self.filename
            with nogil:
                htsfile = hts_open(cfilename, cmode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile

    cpdef Iterator fetch(self, int stid, int maxtid):
        '''fetch reads aligned on chromosome stid~maxtid
        
        Usage: fetch(0,2) to fetch reads aligned on 0,1 (tid)

        Parameters
        ----------
        stid: int
            tid of start chromosome, include.

        maxtid: int
            tid of end chromosome, not include.
        
        Returns
        -------
        Iterator: Iterator
		    An iterator over a collection of reads.
        '''
        return Iterator(self, stid, maxtid)
#
# ---------------------------------------------------------------
#
cdef class IteratorSingle:
    def __cinit__(self):
        self.b = <bam1_t*>calloc(1, sizeof( bam1_t))
        if self.b == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self, BamFile bamfile, int tid):
        self.htsfile = bamfile.htsfile
        self.index = bamfile.index
        self.retval = 0
        with nogil:
            self.iter = sam_itr_queryi(self.index, tid, 0, MAX_POS)

    cdef int cnext(self):
        '''cversion of iterator. retval>=0 if success.'''
        with nogil:
            self.retval = hts_itr_next(self.htsfile.fp.bgzf, self.iter, self.b, self.htsfile)

    def __dealloc__(self):
        bam_destroy1(self.b)
        hts_itr_destroy(self.iter)
#
# ---------------------------------------------------------------
#
cdef class Iterator:
    def __init__(self, BamFile bamfile, int stid, int maxtid):
        self.bamfile = bamfile
        self.htsfile = bamfile.htsfile
        self.index = bamfile.index
        self.tid = -1
        self.maxtid = maxtid

    def nextiter(self):
        # get a new iterator for a chromosome. The file will not be re-opened.
        self.rowiter = IteratorSingle(self.bamfile, self.tid)

    def __iter__(self):
        return self

    def __next__(self):
        '''fetch a record and parse it's CIGAR

        Returns
        -------
        Iterator: Iterator
		    An iterator over a collection of reads.
        '''
        cdef list segl = []
        # Create an initial iterator
        if self.tid == -1:
            self.tid = self.stid
            self.nextiter()

        while 1:
            self.rowiter.cnext()
            # If current iterator is not exhausted, return aligned read
            if self.rowiter.retval > 0:
                segl = parse_cigar(self.rowiter.b)
                return segl

            # Otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.maxtid:
                self.nextiter()
            else:
                raise StopIteration