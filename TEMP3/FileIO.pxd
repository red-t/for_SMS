from .HtslibExternal cimport *
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.errno  cimport errno
from libc.string cimport strerror

cdef class BamFile:
    """BamFile(str filePath, str mode, int numThread=1, BamFile template=None)"""

    cdef char       *filePath
    cdef char       *indexFilePath
    cdef char       *mode
    cdef int        numThread
    cdef htsFile    *htsFile
    cdef hts_idx_t  *index
    cdef sam_hdr_t  *header

    cdef void      openFile(self, BamFile template=*)
    cdef htsFile   *openHtsFile(self) except? NULL
    cdef void      write(self, bam1_t *bamRecord)


cdef class Iterator:
    """Iterator(BamFile bamFile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef bam1_t     *bamRcord
    cdef htsFile    *htsFile
    cdef hts_itr_t  *iter
    cdef int64_t    offset

    cdef int cnext1(self)
    cdef int cnext2(self)
    cdef int cnext3(self, int64_t offset)


cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    ctypedef struct AiList:
        pass


cdef extern from "src/cluster_utils.h" nogil:
    ctypedef packed struct Cluster:
        int         refStart
        int         refEnd
        int         startIndex
        int         endIndex
        float       numSeg
        uint16_t    directionFlag
        uint8_t     cltType
        uint8_t     locationType
        uint8_t     numSegType
        float       entropy
        float       balanceRatio
        float       lowMapQualFrac
        float       dualClipFrac
        float       alnFrac1
        float       alnFrac2
        float       alnFrac4
        float       alnFrac8
        float       alnFrac16
        float       meanMapQual
        float       meanAlnScore
        float       meanQueryMapFrac
        float       meanDivergence
        float       bgDiv
        float       bgDepth
        float       bgReadLen
        float       teAlignedFrac
        int         teTid
        uint8_t     isInBlacklist
        float       probability
    
    ctypedef struct Args:
        int         numThread
        int         tid
        int         minSegLen
        int         maxDistance
        int         numTeTid
        int         *teTidCountTable
        int         minOverhang
        float       bgDiv
        float       bgDepth
        float       bgReadLen
        htsFile     *genomeBamFile
        bam1_t      *firstBamRecord
        bam1_t      *secondBamRecord
        AiList      *repeatAiList
        AiList      *gapAiList
        AiList      *blackAiList
    
    int overhangIsShort(Segment *segment, int minOverhang)
    int isLowQualClt(Cluster *cluster)
    int isGermClt(Cluster *cluster)
    int isSomaClt(Cluster *cluster)
    int getOuputSegIndex(Cluster *cluster, Segment *segArray, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct Segment:
        uint16_t    flag
        uint8_t     mapQual
        int         queryStart
        int         queryEnd
        int         refPosition
        uint8_t     segType
        uint8_t     alnType
        int64_t     fileOffset
        int         alnRefStart
        int         alnRefEnd
        uint8_t     order
        uint8_t     numSeg
        int         overhang
        int         matchLen
        int         readLen
        uint8_t     alnLocationType
        uint8_t     numTeAlignment
        int         sumQueryMapLen
        float       sumAlnScore
        float       sumDivergence
        uint16_t    directionFlag
        int         startIndex
        int         endIndex
        int         teTid

    void getTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize)
    int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIndex, int sourceStart, int sourceEnd)


cdef Args newArgs(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs)
cdef ouputAllSegSeqs(Segment[::1] segArray, BamFile genomeBamFile, Args args)
cdef outputGermCltSeqs(Cluster[::1] cltArray, Segment[::1] segArray, BamFile genomeBamFile, Args args)
cpdef outputSomaCltSeqs(Cluster[::1] cltArray, Segment[::1] segArray, object cmdArgs, int tid)