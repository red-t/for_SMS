from AlignmentFileIO cimport BamFile
from SegmentParser cimport *

cdef dict CLUSTER_DICT

cdef class Cluster:
    """a class representing cluster of candidate support reads"""
    cdef public:
        list segments
        int  start, end
        int  nseg
    
    cpdef add(self, InsertSegment iseg)