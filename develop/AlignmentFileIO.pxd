from htslib_external cimport *
from SegmentParser cimport parse_cigar, InsertSegment
from libc.stdlib cimport malloc, calloc, realloc, free
from pysam.libcutils cimport encode_filename, force_str
from libc.errno  cimport errno
from libc.string cimport strerror


cdef class BamFile:
    """BamFile(str filepath, int32_t nthreads=1, str mode='rb', BamFile template=None)

    A Class for reading `BAM` formatted file:
    1. For reading, you don't need to provide an existing BamFile instance as template,
       and for writing, a template is necessary.
    2. Using `BamFile.fetch(mintid, maxtid, minl)` to traverse all aligned reads and e-
       xtracted "insert segments".
    3. Using `BamFile.write(iseg)` to write single alignment of `iseg` to disk.

    Parameters
    ----------
    filepath: str
        Alias for `filename`.
    nthreads: int32_t
        Number of threads to use for decompressing BAM files. (Default=1)
    mode: str
        opening mode, 'rb' for reading, 'wb' for writing. (Default='rb')
    template: BamFile
        an exisiting BamFile instance, it is necessary for writting. (Default=None)
    """
    # cdef    AlignmentHeader header
    cdef    htsFile *htsfile        # pointer to htsFile structure
    cdef    hts_idx_t *index        # pointer to hts_idx_t structure
    cdef    bam_hdr_t  *hdr         # pointer to bam_hdr_t structure
    cdef    char *filename          # filename as supplied by user
    cdef    char *index_filename    # filename of index
    cdef    char *mode              # opening mode
    cdef    int32_t threads         # number of threads to use

    cdef  void      _open(self, BamFile template=*)
    cdef  htsFile   *_open_htsfile(self) except? NULL
    cpdef dict      fetch(self, uint8_t stid, uint8_t maxtid, uint8_t minl=*)
    cpdef void       write(self, InsertSegment iseg)



cdef class IteratorSingle:
    """IteratorSingle(BamFile bamfile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
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