from .FileIO cimport *

cdef extern from "src/cluster_utils.h" nogil:
    void updateCluster(Cluster *cltArray, Segment *segArray, Args args)
    void intersectBlackList(Cluster *cluster, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct TeAlignment:
        int     segIdx
        int     AlnScore
        int     queryStart
        int     queryEnd
        int     mapLen
        float   divergence
    
    int fillSegArray(bam1_t *bam, Segment *segArray, int64_t fileOffset, int minSegmentLength)
    void updateSegment(Segment *segArray, AiList *repeatAiList, AiList *gapAiList)
    void updateSegByTeArray(Segment *segArray, TeAlignment *teArray, int teIdx)

    #########################
    ### Alignment records ###
    #########################
    void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bam)
    void fillTeArray(bam1_t *bam, TeAlignment *teArray)
    

cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid)
cdef object getHighQualClts(dict allCltData)