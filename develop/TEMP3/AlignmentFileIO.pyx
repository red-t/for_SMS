import os
import numpy as np

#
# ---------------------------------------------------------------
#
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
            sam_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            sam_hdr_destroy(self.hdr)
            self.hdr = NULL

    def __dealloc__(self):
        if self.htsfile:
            sam_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            sam_hdr_destroy(self.hdr)
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
        elif (self.mode == b'wb') or (self.mode == b'wF'):
            # copy header from template
            if template:
                self.hdr = sam_hdr_dup(template.hdr)
            else:
                raise ValueError("not enough information to construct header. Please provide template")
            
            # write header to htsfile
            if self.mode == b'wb':
                with nogil:
                    sam_hdr_write(self.htsfile, self.hdr)
    
    cdef htsFile *_open_htsfile(self) except? NULL:
        '''open file in 'rb/wb' mode, return htsFile object if success.'''
        cdef int32_t threads = self.threads
        cdef htsFile *htsfile

        if isinstance(self.filename, bytes):
            with nogil:
                htsfile = sam_open(self.filename, self.mode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile

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
            raise IOError("sam_write1 failed with error code {}".format(ret))
#
# ---------------------------------------------------------------
#
cdef class Iterator:
    def __cinit__(self):
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))
        if self.b == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self,
                 BamFile bamfile,
                 int tid):

        self.bamfile = bamfile
        self.htsfile = bamfile.htsfile
        self.index   = bamfile.index
        self.tid     = tid
        with nogil:
            self.iter = sam_itr_queryi(self.index,
                                       tid,
                                       0,
                                       MAX_POS)
    
    def __dealloc__(self):
        bam_destroy1(self.b)
        sam_itr_destroy(self.iter)
    
    def __iter__(self):
        return self
    
    cdef int cnext_offt(self, int64_t offset):
        '''read a alignment record with specified offset'''
        cdef int retval

        with nogil:
            bgzf_seek(self.htsfile.fp.bgzf, offset, SEEK_SET)
            retval = bam_read1(self.htsfile.fp.bgzf,
                               self.b)
        
        return retval

    cdef int cnext(self):
        '''cversion of iterator. retval>=0 if success.'''
        cdef int        retval
        cdef int64_t    offset

        with nogil:
            offset = bgzf_tell(self.htsfile.fp.bgzf)
            retval = hts_itr_next(self.htsfile.fp.bgzf,
                                  self.iter,
                                  self.b,
                                  self.htsfile)
        self.offset = offset
        return retval