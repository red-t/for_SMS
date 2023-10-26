#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H
#include <stdint.h>
#include <math.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "seg_utils.h"
#include "AIList.h"
//-------------------------------------------------------------------------------------

/******************
 *** Structures ***
 ******************/
/*! @typedef
 @abstract Structure cluster merged from segments.
 @field  refStart           cluster start on reference sequence (0-based, included)
 @field  refEnd             cluster end on reference sequence (0-based, not-included)
 @field  startIndex         start index in segments array (0-based, include)
 @field  endIndex           start index in segments array (0-based, not-include)
 @field  numSeg             number of segments in the cluster (normalized by bg depth)
 @field  directionFlag      bitwise flag representing cluster direction
                                1: forward
                                2: reverse
 @field  cltType            cluster type
                                0: germline (multiple support reads)
                                1: somatic (1 support read & 1 alignment)
                                2: somatic (1 support read & 2 alignments)
 @field  locationType       bitwise flag representing cluster location
                                1: inside normal region
                                2: at repeat/gap boundary
                                4: inside repeat/gap
 @field  numSegType         number of different segment types
 @field  entropy            entropy based on fraction of different type segments
 @field  balanceRatio       balance ratio based on number of left- & right-clip segments
 @field  lowMapQualFrac     fraction of segments with low mapQual (<5)
 @field  dualClipFrac       fraction of "dual-clip" alignments
 @field  alnFrac1           fraction of segments with alnLocationType=1
 @field  alnFrac2           fraction of segments with alnLocationType=2
 @field  alnFrac4           fraction of segments with alnLocationType=4
 @field  alnFrac8           fraction of segments with alnLocationType=8
 @field  alnFrac16          fraction of segments with alnLocationType=16
 @field  meanMapQual        mean mapQual of cluster
 @field  meanAlnScore       mean per-base alignment score (based on teAlignments)
 @field  meanQueryMapFrac   mean query mapped fraction (based on teAlignments)
 @field  meanDivergence     mean per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
 @field  bgDiv              background divergence (for normalization)
 @field  bgDepth            background depth (for normalization)
 @field  bgReadLen          background read length
 @field  teTid              majority TE-tid of cluster
 */
typedef struct {
    int         refStart;
    int         refEnd;
    int         startIndex;
    int         endIndex;
    float_t     numSeg;
    uint16_t    directionFlag;
    uint8_t     cltType;
    uint8_t     locationType;
    uint8_t     numSegType;
    float_t     entropy;
    float_t     balanceRatio;
    float_t     lowMapQualFrac;
    float_t     dualClipFrac;
    float_t     alnFrac1;
    float_t     alnFrac2;
    float_t     alnFrac4;
    float_t     alnFrac8;
    float_t     alnFrac16;
    float_t     meanMapQual;
    float_t     meanAlnScore;
    float_t     meanQueryMapFrac;
    float_t     meanDivergence;
    float_t     bgDiv;
    float_t     bgDepth;
    float_t     bgReadLen;
    float_t     teAlignedFrac;
    int         teTid;
} __attribute__((packed)) Cluster;

typedef struct {
    int         numThread;
    int         tid;
    int         minSegLen;
    int         maxDistance;
    int         numTeTid;
    int         *teTidCountTable;
    int         minOverhang;
    float       bgDiv;
    float       bgDepth;
    float       bgReadLen;
    htsFile     *genomeBamFile;
    bam1_t      *firstBamRecord;
    bam1_t      *secondBamRecord;
    AiList      *repeatAiList;
    AiList      *gapAiList;
} Args;


/**********************
 *** Update Cluster ***
 **********************/
Args initArgs(int numThread, int tid, int minSegLen, int maxDistance, int minOverhang, float bgDiv, float bgDepth, float bgReadLen);
void updateCluster(Cluster *cltArray, Segment *segArray, Args args);


/************************
 *** Define Candidate ***
 ************************/
#define overhangIsShort(segment, minOverhang) ((segment)->overhang < (minOverhang))
#define nameIsSame(record1, record2) (strcmp(bam_get_qname((record1)), bam_get_qname((record2))) == 0)
#define isValidCandidate(cluster) ((cluster)->teAlignedFrac >= 0.8)

void setTeAlignedFrac(Cluster *cluster, Segment *segArray, Args args);
void setCltType(Cluster *cluster, Segment *segArray, Args args);


/**********************************
 *** Update Cluster By Segments ***
 **********************************/
#define isDualClip(segment) (((segment)->alnType & DUAL_CLIP) == 5)
#define noTeAlignment(segment) ((segment)->numTeAlignment == 0)
#define directionIsConsistent(record) (((record)->directionFlag & 255) > ((record)->directionFlag >> 8))
#define directionIsInconsistent(record) (((record)->directionFlag & 255) < ((record)->directionFlag >> 8))

void updateBySegArray(Cluster *cluster, Segment *segArray, Args args);
void countValuesFromSeg(Cluster *cluster, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight);
static inline void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment);
static inline void countAlnFracs(Cluster *cluster, Segment *segment);
static inline void setEntropy(Cluster *cluster, int numLeft, int numMiddle, int numRight);
float_t getEntropy(int numLeft, int numMiddle, int numRight, int numSeg);
static inline void setBalanceRatio(Cluster *cluster, int numLeft, int numRight);
static inline void setNumSegType(Cluster *cluster);


/*********************************
 *** Update By Backbg Info ***
 *********************************/
void setCltLocationType(Cluster *cluster, Args args);
static inline void divideValuesByNumSeg(Cluster *cluster);
static inline void divideValuesByBackbg(Cluster *cluster, Args args);
static inline void setDirection(Cluster *cluster);
static inline void setBackbgInfo(Cluster *cluster, Args args);

#endif // CLUSTER_UTILS_H