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

    # Add a interval into AIList
    # @param s   start of the interval
    # @param e   end of the interval
    void ailist_add(ailist_t *ail, uint32_t s, uint32_t e)

    # Add intervals from BED-like file
    # @param fn   filename of the BED-like file
    # @param chr  name of the specified chromosome, only intervals
    #             on this chromosome will be added
    void readBED(ailist_t *ail, const char* fn, char* chr)

    # Construct ailist: decomposition and augmentation
    # @param cLen minimum coverage length, default=20
    void ailist_construct(ailist_t *ail, int cLen)

    # Compute the minimum distance between query position to the nearest interval
    # @param rpos  reference position used as query
    # @param flank flank size used for extending query
    # @param d     minimum distance, will be computed when calling this function
    # @return      number of intervals overlapped with [rpos-flank, rpos+flank)
    uint32_t query_dist_p(ailist_t *ail, uint32_t rpos, uint32_t flank, uint32_t *d)

    # Compute the minimum distance between query interval to the nearest interval
    # @param st    start position of query interval
    # @param ed    end position of query interval
    # @param flank flank size used for extending query interval
    # @param d     minimum distance, will be computed when calling this function
    # @return      number of intervals overlapped with [st-flank, ed+flank)
    uint32_t query_dist_c(ailist_t *ail, uint32_t st, uint32_t ed, uint32_t flank, uint32_t *d)

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