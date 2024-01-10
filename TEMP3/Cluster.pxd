from .FileIO cimport *
from .HtslibExternal cimport *

cdef extern from "src/cluster_utils.h" nogil:
    int isValidCandidate(Cluster *cluster)
    void updateCluster(Cluster *cltArray, Segment *segArray, Args args)
    void intersectBlackList(Cluster *cluster, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct TeAlignment:
        int     segIndex
        int     AlnScore
        int     queryStart
        int     queryEnd
        int     mapLen
        float   divergence
        int16_t flag
        int     teTid
    
    int fillSegmentArray(bam1_t *bamRecord, Segment *segmentArray, int64_t fileOffset, int minSegmentLength)
    void updateSegment(Segment *segArray, AiList *repeatAiList, AiList *gapAiList)
    void updateSegByTeArray(Segment *segArray, TeAlignment *teArray, int teIndex)
    void countTeTids(Segment *segArray, TeAlignment *teArray, int *teTidCountTable, int numTeTid)

    #########################
    ### Alignment records ###
    #########################
    int bamIsInvalid(bam1_t *bamRecord)
    void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bamRecord)
    void fillTeArray(bam1_t *bamRecord, TeAlignment *teArray)


cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    AiList *initAiList()
    void destroyAiList(AiList *ail)
    void readBED(AiList *ail, const char* bedFileName, const char* targetChrom)
    void constructAiList(AiList *ail, int minCoverageLen)

cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid)
cdef object getHighQualClts(dict chromCltData)