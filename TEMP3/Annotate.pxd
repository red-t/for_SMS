from .FileIO cimport *

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Anno:
        int     idx
        int     cltTid
        int     cltIdx
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
    void annoTsd(Cluster *cluster)
    void outPutAnno(Anno *annoArray, int numAnno, const char *teFn, const char *outFn)

cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs)