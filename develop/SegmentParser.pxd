from htslib_external cimport *
from pysam.libcutils cimport force_bytes
from cpython cimport PyBytes_FromStringAndSize

cdef int parse_cigar(bam1_t *src,
                     int[:, ::1] segs,
                     int N,
                     int offset,
                     int minl=*)

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
                     uint8_t minl=*)