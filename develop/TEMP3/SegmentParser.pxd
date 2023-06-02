from .htslib_external cimport *
from cpython cimport PyBytes_FromStringAndSize

cdef extern from "src/seg_utils.h" nogil:
    ######################
    ### CIGAR resolver ###
    ######################

    # Get whether the CIGAR operation is clip/insert.
    # @param  op CIGAR operation
    # @return    1 if the CIGAR operation is clip/insert, 0 if not
    int is_clip_or_insert(uint32_t op)

    # Get whether the CIGAR operation is match/equal/diff.
    # @param  op CIGAR operation
    # @return    1 if the CIGAR operation is match/equal/diff, 0 if not
    int is_match(uint32_t op)

    # Get whether the CIGAR operation is del/skip.
    # @param  op CIGAR operation
    # @return    1 if the CIGAR operation is del/skip, 0 if not
    int is_del_or_skip(uint32_t op)


    #######################
    ### Segment records ###
    #######################

    uint8_t LEFT_CLIP
    uint8_t MID_INSERT
    uint8_t RIGHT_CLIP
    uint8_t DUAL_CLIP

    ctypedef packed struct seg_dtype_struct:
        uint16_t    flag
        uint8_t     mapq
        int32_t     qst
        int32_t     qed
        int32_t     rpos
        int32_t     lqseq
        uint8_t     sflag
        uint8_t     rflag
        int64_t     offset
        int32_t     refst
        int32_t     refed
        uint8_t     ith
        uint8_t     nseg
        int32_t     overhang
        int32_t     nmatch
        uint8_t     loc_flag

    # # Get whether the alignment is dual-clip.
    # # @param rflag  rflag of the segment, representing type of the corresponding alignment
    # # @return       1 if the alignment is dual-clip, 0 if not
    # int aln_is_dclip(uint8_t rflag)

    # # Get whether the query is secondary alignment by flag.
    # # @param flag  bitwise flag of the query alignment
    # # @return      1 if query is secondary, 0 if not
    # int aln_is_second(uint16_t flag)

    # # Get whether the query is on the reverse strand by flag.
    # # @param flag  bitwise flag of the query alignment
    # # @return      1 if query is on the reverse strand, 0 if not
    # int aln_is_rev(uint16_t flag)

    # # Update overhang for mid-insert type segment
    # # @param overhang  original overhang of the segment
    # # @param nmatch    nummber of match bases of the corresponding alignment
    # # @return          Updated overhang, <= original overhang
    # int32_t update_overhang(int32_t overhang, int32_t nmatch)


cdef int parse_cigar(bam1_t *src,
                     seg_dtype_struct[::1] segs,
                     int32_t N,
                     int64_t offset,
                     int32_t minl=*)
#
# ---------------------------------------------------------------
#
#####################
### Previous code ###
#####################
# cdef class InsertSegment:
#     """A class for insert/clip segment"""
#     cdef bam1_t *_delegate
#     cdef public:
#         int32_t qstart, qend        # position of the segment on query sequence, [qstart, qend) region of query sequence
#         int32_t rpos                # position of the segment on reference sequence, rend=rpos, rstart=rpos-1
#         int16_t stype               # left-clip: 0x1, mid-insert: 0x2, right-clip: 0x4
#         int32_t ref_end             # reference position of the alignment, ref_start=_delegate.core.pos
#         int32_t q_start, q_end      # query position of the alignment
#         int32_t overhang            # minimum length of the anchor part
#         int32_t ldist, rdist        # left/right distance to the neighbor segment
    
#     cpdef int64_t get_tag_i(self, str tag)
#     cpdef double  get_tag_f(self, str tag)
#     cpdef str     get_seq(self, int start=*, int end=*)
#     cpdef tuple   trim(self, int tsize)
# #
# # ---------------------------------------------------------------
# #
# cdef InsertSegment makeInsertSegment(bam1_t *src,
#                                      int32_t qstart,
#                                      int32_t qend,
#                                      int32_t rpos,
#                                      int16_t stype)
# #
# # ---------------------------------------------------------------
# #
# cdef int parse_cigar1(bam1_t *src,
#                      list tmp_segl,
#                      int32_t minl=*)