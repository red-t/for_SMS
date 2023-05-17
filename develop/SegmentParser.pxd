from htslib_external cimport *
from pysam.libcutils cimport force_bytes
from cpython cimport PyBytes_FromStringAndSize

cdef packed struct seg_dtype_struct:
    uint16_t    flag
    uint8_t     mapq
    int32_t     qst
    int32_t     qed
    int32_t     rpos
    uint8_t     sflag
    uint8_t     rflag
    int64_t     offset
    int32_t     refst
    int32_t     refed
    uint8_t     ith
    uint8_t     nseg
    int32_t     overhang
    int32_t     nmatch
    uint8_t     loc_flag1
    uint8_t     loc_flag2

cdef int parse_cigar(bam1_t *src,
                     seg_dtype_struct[::1] segs,
                     int32_t N,
                     int64_t offset,
                     int32_t minl=*)

cdef class InsertSegment:
    """A class for insert/clip segment"""
    cdef bam1_t *_delegate
    cdef public:
        int32_t qstart, qend        # position of the segment on query sequence, [qstart, qend) region of query sequence
        int32_t rpos                # position of the segment on reference sequence, rend=rpos, rstart=rpos-1
        int16_t stype               # left-clip: 0x1, mid-insert: 0x2, right-clip: 0x4
        int32_t ref_end             # reference position of the alignment, ref_start=_delegate.core.pos
        int32_t q_start, q_end      # query position of the alignment
        int32_t overhang            # minimum length of the anchor part
        int32_t ldist, rdist        # left/right distance to the neighbor segment
    
    cpdef int64_t get_tag_i(self, str tag)
    cpdef double  get_tag_f(self, str tag)
    cpdef str     get_seq(self, int start=*, int end=*)
    cpdef tuple   trim(self, int tsize)



cdef InsertSegment makeInsertSegment(bam1_t *src,
                                     int32_t qstart,
                                     int32_t qend,
                                     int32_t rpos,
                                     int16_t stype)


cdef int parse_cigar1(bam1_t *src,
                     list tmp_segl,
                     int32_t minl=*)