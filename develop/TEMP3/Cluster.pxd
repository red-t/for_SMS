from .AlignmentFileIO cimport BamFile
from .SegmentParser cimport seg_dtype_struct
from .htslib_external cimport *


cdef packed struct cluster_dtype_struct:
    int32_t     st
    int32_t     st_idx
    int32_t     ed
    int32_t     ed_idx
    int16_t     nseg
    uint8_t     strand
    uint8_t     cloc_flag
    uint8_t     n1
    int32_t     l1
    uint8_t     n2
    int32_t     l2
    uint8_t     ntype
    float       entropy
    float       bratio
    int16_t     nL
    int16_t     nM
    int16_t     nR
    float       lmq_freq
    float       avg_mapq
#
# ---------------------------------------------------------------
#
##############
### AIList ###
###############

cdef extern from "src/AIList.h" nogil:
    ctypedef struct ctg_t:
        pass

    ctypedef struct ailist_t:
        ctg_t *ctg

    # Initialize ailist_t
    ailist_t *ailist_init()

    # Free ailist data
    void ailist_destroy(ailist_t *ail)

    # Add intervals from BED-like file
    # @param fn   filename of the BED-like file
    # @param chr  name of the specified chromosome, only intervals
    #             on this chromosome will be added
    void readBED(ailist_t *ail, const char* fn, char* chr)

    # Construct ailist: decomposition and augmentation
    # @param cLen minimum coverage length, default=20
    void ailist_construct(ailist_t *ail, int cLen)

    # Compute alignment location flag
    # @param rep_ail	AIList of repeats
    # @param gap_ail	AIList of gaps
    # @param segs		address to the array of segments
    void aln_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, seg_dtype_struct *)

#
# ---------------------------------------------------------------
#
cdef object merge_segments(seg_dtype_struct[::1] segs,
                                            int maxdist)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str fpath,
                         int threads,
                         int tid,
                         int minl,
                         int maxdist)