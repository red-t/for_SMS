from .htslib_external cimport *
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.errno  cimport errno
from libc.string cimport strerror

cdef class BamFile:
    """BamFile(str filePath, str mode, int numThread=1, BamFile template=None)"""

    cdef char       *filePath
    cdef char       *indexFilePath
    cdef char       *mode
    cdef int        numThread
    cdef htsFile    *htsFile
    cdef hts_idx_t  *index
    cdef sam_hdr_t  *header

    cdef void      openFile(self, BamFile template=*)
    cdef htsFile   *openHtsFile(self) except? NULL
    cdef void      write(self, bam1_t *bamRecord)


cdef class Iterator:
    """Iterator(BamFile bamFile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef bam1_t     *bamRcord
    cdef htsFile    *htsFile
    cdef hts_itr_t  *iter
    cdef int64_t    offset

    cdef int cnext1(self)
    cdef int cnext2(self)
    cdef int cnext3(self, int64_t offset)