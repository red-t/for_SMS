from htslib_external cimport *
from SegmentParser cimport parse_cigar
from libc.stdlib cimport malloc, calloc, realloc, free
from pysam.libcutils cimport encode_filename, force_bytes, force_str
from libc.errno  cimport errno
from libc.string cimport strerror


cdef class BamFile:
    """BamFile(str filepath, int nthreads=1)

    A Class for reading `BAM` formatted file. Using `BamFile.fetch(mintid, maxtid)`
    to traverse all aligned reads and parse their CIGAR.

    Parameters
    ----------
    filepath : str
        Alias for `filename`.
    
    threads: int
        Number of threads to use for decompressing BAM files. (Default=1)
    """
    # cdef    AlignmentHeader header
    cdef    htsFile *htsfile        # pointer to htsFile structure
    cdef    hts_idx_t *index        # pointer to hts_idx_t structure
    cdef    char *filename          # filename as supplied by user
    cdef    char *index_filename    # filename of index
    cdef    int32_t threads         # number of threads to use

    cdef _open(self)
    cdef htsFile *_open_htsfile(self) except? NULL
    cpdef dict fetch(self, uint8_t stid, uint8_t maxtid, uint8_t minl=*)



cdef class IteratorSingle:
    """IteratorSingle(BamFile bamfile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef int retval
    cdef bam1_t *b          # pointer to a record in BAM file, change when call `cnext`
    cdef htsFile *htsfile   # pointer to htsFile structure
    cdef hts_idx_t *index   # pointer to hts_idx_t structure
    cdef hts_itr_t *iter    # pointer to hts_itr_t structure, iterator from htslib

    cdef int cnext(self)



cdef class Iterator:
    """Iterator(BamFile bamfile, int stid, int maxtid)

    A class for iterating over mapped reads in all specified chromosomes
    """
    cdef BamFile bamfile
    cdef htsFile *htsfile
    cdef hts_idx_t *index
    cdef IteratorSingle rowiter
    cdef uint8_t stid, maxtid, minl
    cdef int8_t tid