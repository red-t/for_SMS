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
                         segs)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str     fpath,
                         int32_t threads,
                         uint8_t tid,
                         uint8_t minl,
                         uint8_t maxdist)