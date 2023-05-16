import os

cdef class BamFile:
    def __cinit__(self,
                  str     filepath,
                  int32_t nthreads = 1,
                  str     mode     = 'rb',
                  BamFile template = None):
        cdef bytes bfilepath    = os.fsencode(filepath)
        cdef bytes bidx         = os.fsencode(filepath+".bai")
        cdef bytes bmode        = mode.encode(TEXT_ENCODING, ERROR_HANDLER)

        self.threads    = nthreads
        self.mode       = bmode
        self.filename   = bfilepath
        self.index_filename = bidx
        self._open(template)
        
    def close(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            bam_hdr_destroy(self.hdr)
            self.hdr = NULL

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            bam_hdr_destroy(self.hdr)
            self.hdr = NULL

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    cdef void _open(self, BamFile template=None):
        '''open BAM file for reading/writting'''
        # open file as htsFile 
        self.htsfile = self._open_htsfile()
        if self.htsfile == NULL:
            if errno:
                raise IOError(errno, "could not open alignment file `{}`: {}".format(
                    self.filename.decode(TEXT_ENCODING, ERROR_HANDLER),
                    strerror(errno)))
            else:
                raise ValueError("could not open alignment file `{}`".format(
                    self.filename.decode(TEXT_ENCODING, ERROR_HANDLER)))

        # for reading
        if self.mode == b'rb':
            # bam files require a valid header
            with nogil:
                self.hdr = sam_hdr_read(self.htsfile)
            if self.hdr == NULL:
                raise ValueError("file does not have a valid header, is it BAM format?")
            
            # open BAM index file to enable random access
            with nogil:
                self.index = sam_index_load2(self.htsfile, self.filename, self.index_filename)
            if not self.index:
                if errno:
                    raise IOError(errno, strerror(errno))
                else:
                    raise IOError('unable to open index file `{}`'.format(
                        self.index_filename.decode(TEXT_ENCODING, ERROR_HANDLER)))
        # for writing
        elif self.mode == b'wb':
            # copy header from template
            if template:
                self.hdr = bam_hdr_dup(template.hdr)
            else:
                raise ValueError("not enough information to construct header. Please provide template")
            
            # write header to htsfile
            with nogil:
                sam_hdr_write(self.htsfile, self.hdr)
    
    cdef htsFile *_open_htsfile(self) except? NULL:
        '''open file in 'rb/wb' mode, return htsFile object if success.'''
        cdef int32_t threads = self.threads
        cdef htsFile *htsfile

        if isinstance(self.filename, bytes):
            with nogil:
                htsfile = hts_open(self.filename, self.mode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile

    cpdef dict fetch(self,
                     BamFile wbf,
                     uint8_t stid,
                     uint8_t maxtid,
                     uint8_t minl=50):
        '''fetch reads aligned on chromosome stid ~ maxtid-1
        
        Usage: fetch(0, 2) to fetch reads aligned on 0,1 (tid)

        Parameters
        ----------
        wbf: BamFile
            BamFile object opened for writting.

        stid: uint8_t
            tid of start chromosome, include.

        maxtid: uint8_t
            tid of end chromosome, not include.

        minl: uint8_t
            minimum length of the segment. segment with cigar_length < minl
            will not be used to create a new InsertSegment.
        
        Returns
        -------
        SEG_DICT: dict
		    dict of list, tid -> list_of_InsertSegment
        '''
        cdef Iterator ite
        cdef int i

        for i in range(stid, maxtid):
            SEG_DICT[i] = []
        
        ite = Iterator(self, wbf, stid, maxtid, minl)
        for i in ite:
            continue

        return SEG_DICT
    
    cdef void write(self, bam1_t *src):
        '''write a single alignment to disk.

        Parameters
        ----------
        src: bam1_t*
            bam1_t type pointer, which point to the source memory address
            of the alignment loaded into memory.
        '''
        cdef int ret

        with nogil:
            ret = sam_write1(self.htsfile, self.hdr, src)
        if ret < 0:
            raise IOError(
            "sam_write1 failed with error code {}".format(ret))
#
# ---------------------------------------------------------------
#
cdef class IteratorSingle:
    def __cinit__(self):
        self.b = <bam1_t*>calloc(1, sizeof( bam1_t))
        if self.b == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self,
                 BamFile bamfile,
                 int tid):

        self.htsfile = bamfile.htsfile
        self.index   = bamfile.index
        with nogil:
            self.iter = sam_itr_queryi(self.index,
                                       tid,
                                       0,
                                       MAX_POS)

    cdef int cnext(self):
        '''cversion of iterator. retval>=0 if success.'''
        cdef int retval
        with nogil:
            retval = hts_itr_next(self.htsfile.fp.bgzf,
                                  self.iter,
                                  self.b,
                                  self.htsfile)
            return retval

    def __dealloc__(self):
        bam_destroy1(self.b)
        hts_itr_destroy(self.iter)
#
# ---------------------------------------------------------------
#
cdef class Iterator:
    def __init__(self,
                 BamFile bamfile,
                 BamFile wbf,
                 uint8_t stid,
                 uint8_t maxtid,
                 uint8_t minl):
                 
        self.bamfile = bamfile
        self.wbf     = wbf
        self.htsfile = bamfile.htsfile
        self.index   = bamfile.index
        self.tid     = -1
        self.stid    = stid
        self.maxtid  = maxtid
        self.minl    = minl

    def nextiter(self):
        # get a new iterator for a chromosome. The file will not be re-opened.
        self.rowiter = IteratorSingle(self.bamfile, self.tid)

    def __iter__(self):
        return self

    def __next__(self):
        '''fetch a record and parse it's CIGAR, store results in SEG_DICT'''
        cdef int    n, l, retval
        cdef list   tmp_segl = []
        
        # Create an initial iterator
        if self.tid == -1:
            self.tid = self.stid
            self.nextiter()

        while 1:
            retval = self.rowiter.cnext()
            # If current iterator is not exhausted, return aligned read
            if retval > 0:
                n = self.rowiter.b.core.n_cigar
                l = self.rowiter.b.core.l_qseq
                if n==0 or l==0:
                    continue
                retval = parse_cigar(self.rowiter.b,
                                     tmp_segl,
                                     self.minl)
                if retval > 0:
                    self.wbf.write(self.rowiter.b)
                continue

            SEG_DICT[self.tid] = tmp_segl
            tmp_segl = []
            # Otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.maxtid:
                self.nextiter()
            else:
                raise StopIteration