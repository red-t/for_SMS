from AlignmentFileIO cimport BamFile
from SegmentParser cimport *

cdef dict CLUSTER_DICT
#
# ---------------------------------------------------------------
#
cdef class Cluster:
    """a class representing cluster of candidate support reads"""
    cdef public:
        list segments
        int  start, end
        int  nseg
    
    cpdef add(self, InsertSegment iseg)
#
# ---------------------------------------------------------------
#
cdef void merge_segments(list    segments,
                         int     tid,
                         int     maxdist,
                         object segs,
                         object clts)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str     fpath,
                         int32_t threads,
                         uint8_t tid,
                         uint8_t minl,
                         uint8_t maxdist)



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
cdef object merge_segments1(seg_dtype_struct[::1] segs,
                                            int maxdist)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster1(str fpath,
                         int threads,
                         int tid,
                         int minl,
                         int maxdist)