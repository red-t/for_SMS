from .FileIO cimport *

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Anno:
        int     idx
        int     refStart
        int     refEnd
        int     tid
        uint8_t strand
        int     queryStart
        int     queryEnd

cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs)