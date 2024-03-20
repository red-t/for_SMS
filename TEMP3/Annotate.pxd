from .FileIO cimport *

cdef extern from "src/io_utils.h":
    ###################
    ### Cluster I/O ###
    ###################
    void outputClt(Cluster *cltArray, int startIdx, int endIdx, const char *refFn, const char *teFn)

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Anno:
        int         idx
        int         cltTid
        int         cltIdx
        int         queryStart
        int         queryEnd
        uint8_t     strand
        int         tid
        int         refStart
        int         refEnd
        uint32_t    flag
    
    ###################################
    ### Annotate Insertion sequence ###
    ###################################
    int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx)
    void annoTsd(Cluster *cluster, Anno *annoArray, int numAnno)
    void checkGap(Cluster *clt, Anno *annoArray, int numAnno)

    ######################
    ### Annotation I/O ###
    ######################
    void outputAnno(Anno *annoArray, int numAnno, int startIdx, const char *teFn)

cpdef annotateCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs)