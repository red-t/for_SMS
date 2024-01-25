from .FileIO cimport *

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Anno:
        int     idx
        int     queryStart
        int     queryEnd
        uint8_t strand
        int     tid
        int     refStart
        int     refEnd
    
    ###################################
    ### Annotate Insertion sequence ###
    ###################################
    int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx)

cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs)