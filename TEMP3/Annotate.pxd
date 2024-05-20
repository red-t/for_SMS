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
    int fillAnnoArray(Cluster *clt, Anno *annoArray, int idx)
    void annoTsd(Cluster *clt, Anno *annoArray, int numAnno)
    void setInsStruc(Cluster *clt, Anno *annoArray, int numAnno, uint32_t *classArray, int *sizeArray)

    ######################
    ### Annotation I/O ###
    ######################
    void outputAnno(Anno *annoArray, int numAnno, int startIdx, const char *teFn)

cdef extern from "src/post_filter.h":
    ###################
    ### Cluster I/O ###
    ###################
    void postFilter(Cluster *clt)

cpdef annotateCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs)