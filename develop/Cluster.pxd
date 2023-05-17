from AlignmentFileIO cimport BamFile
from SegmentParser cimport *


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