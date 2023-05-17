from htslib_external cimport *
from SegmentParser cimport parse_cigar, parse_cigar1, seg_dtype_struct
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.errno  cimport errno
from libc.string cimport strerror

cdef class BamFile:
    """BamFile(str filepath, int nthreads=1, str mode='rb', BamFile template=None)

    A Class for reading `BAM` formatted file:
    1. For reading, you don't need to provide an existing BamFile instance as template,
       and for writing, a template is necessary.
    2. Using `BamFile.extract_seg(tid, minl)` to traverse all aligned reads and extrac-
       ted "insert segments".
    3. Using `BamFile.write()` to write single alignment to disk.

    Parameters
    ----------
    filepath: str
        Alias for `filename`.

    nthreads: int
        Number of threads to use for decompressing BAM files. (Default=1)

    mode: str
        opening mode, 'rb' for reading, 'wb' for writing. (Default='rb')

    template: BamFile
        an exisiting BamFile instance, it is necessary for writting. (Default=None)
    """
    cdef    htsFile   *htsfile          # pointer to htsFile structure
    cdef    hts_idx_t *index            # pointer to hts_idx_t structure
    cdef    bam_hdr_t *hdr              # pointer to bam_hdr_t structure
    cdef    char      *filename         # filename as supplied by user
    cdef    char      *index_filename   # filename of index
    cdef    char      *mode             # opening mode
    cdef    int       threads           # number of threads to use

    cdef    void      _open(self, BamFile template=*)
    cdef    htsFile   *_open_htsfile(self) except? NULL
    cpdef   dict      extract_seg(self,
                                  BamFile wbf,
                                  int tid,
                                  int minl=*)
    cdef    void       write(self, bam1_t *src)


cdef class Iterator:
    """IteratorSingle(BamFile bamfile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef BamFile    bamfile
    cdef BamFile    wbf
    cdef bam1_t    *b       # pointer to a record in BAM file, change when call `cnext`
    cdef htsFile   *htsfile # pointer to htsFile structure
    cdef hts_idx_t *index   # pointer to hts_idx_t structure
    cdef hts_itr_t *iter    # pointer to hts_itr_t structure, iterator from htslib
    cdef object     segs
    cdef int64_t    offset
    cdef int32_t    minl
    cdef int        tid

    cdef int cnext(self)