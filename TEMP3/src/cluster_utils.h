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
 @field  tid                target id of corresponding chromosome
 @field  refStart           cluster start on reference sequence (0-based, included)
 @field  refEnd             cluster end on reference sequence (0-based, not-included)
 @field  idx                cluster index in clusters array (0-based)
 @field  startIndex         start index in segments array (0-based, include)
 @field  endIndex           end index in segments array (0-based, not-include)
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
 @field  teAlignedFrac      fraction of TE-aligned segments
 @field  teTid              majority TE-tid of cluster
 @field  isInBlacklist      whether cluster intersects with blacklist
 @field  probability        the probability of the cluster to be a positive insertion
 */
typedef struct {
    int         tid;
    int         refStart;
    int         refEnd;
    int         idx;
    int         startIndex;
    int         endIndex;
    float       numSeg;
    uint16_t    directionFlag;
    uint8_t     cltType;
    uint8_t     locationType;
    uint8_t     numSegType;
    float       entropy;
    float       balanceRatio;
    float       lowMapQualFrac;
    float       dualClipFrac;
    float       alnFrac1;
    float       alnFrac2;
    float       alnFrac4;
    float       alnFrac8;
    float       alnFrac16;
    float       meanMapQual;
    float       meanAlnScore;
    float       meanQueryMapFrac;
    float       meanDivergence;
    float       bgDiv;
    float       bgDepth;
    float       bgReadLen;
    float       teAlignedFrac;
    int         teTid;
    uint8_t     isInBlacklist;
    float       probability;
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
    AiList      *blackAiList;
} Args;


/**********************
 *** Update Cluster ***
 **********************/
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

static inline void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment)
{
    switch (segment->segType)
    {
        case LEFT_CLIP:
            *numLeft += 1; break;
        case RIGHT_CLIP:
            *numRight += 1; break;
        case MID_INSERT:
            *numMiddle += 1; break;
        default:
            break;
    }
}

static inline void countAlnFracs(Cluster *cluster, Segment *segment)
{
    switch (segment->alnLocationType)
    {
        case 1:
            cluster->alnFrac1 += 1; break;
        case 2:
            cluster->alnFrac2 += 1; break;
        case 4:
            cluster->alnFrac4 += 1; break;
        case 8:
            cluster->alnFrac8 += 1; break;
        case 16:
            cluster->alnFrac16 += 1; break;
        default:
            break;
    }
}

float getEntropy(int numLeft, int numMiddle, int numRight, int numSeg);

static inline void setEntropy(Cluster *cluster, int numLeft, int numMiddle, int numRight)
{
    cluster->entropy = getEntropy(numLeft, numMiddle, numRight, cluster->numSeg);
}

static inline void setBalanceRatio(Cluster *cluster, int numLeft, int numRight)
{
    cluster->balanceRatio = (MIN(numLeft, numRight) + 0.01) / (MAX(numLeft, numRight) + 0.01);
}

static inline void setNumSegType(Cluster *cluster)
{
    cluster->numSegType = (cluster->numSegType & 1) + ((cluster->numSegType & 2) >> 1) + ((cluster->numSegType & 4) >> 2);
}


/*********************************
 *** Update By Background Info ***
 *********************************/
void setCltLocationType(Cluster *cluster, Args args);

static inline void divideValuesByNumSeg(Cluster *cluster)
{
    if (cluster->alnFrac1 > 0) cluster->alnFrac1 = cluster->alnFrac1 / cluster->numSeg;
    if (cluster->alnFrac2 > 0) cluster->alnFrac2 = cluster->alnFrac2 / cluster->numSeg;
    if (cluster->alnFrac4 > 0) cluster->alnFrac4 = cluster->alnFrac4 / cluster->numSeg;
    if (cluster->alnFrac8 > 0) cluster->alnFrac8 = cluster->alnFrac8 / cluster->numSeg;
    if (cluster->alnFrac16 > 0) cluster->alnFrac16 = cluster->alnFrac16 / cluster->numSeg;

    cluster->dualClipFrac = cluster->dualClipFrac / cluster->numSeg;
    cluster->lowMapQualFrac = cluster->lowMapQualFrac / cluster->numSeg;
    cluster->meanQueryMapFrac = cluster->meanQueryMapFrac / cluster->numSeg;
    cluster->meanDivergence = cluster->meanDivergence / cluster->numSeg;
    cluster->meanMapQual = cluster->meanMapQual / cluster->numSeg;
    cluster->meanAlnScore = cluster->meanAlnScore / cluster->numSeg;
}

static inline void divideValuesByBackbg(Cluster *cluster, Args args)
{
    cluster->meanDivergence = cluster->meanDivergence / args.bgDiv;
    cluster->numSeg = cluster->numSeg / args.bgDepth;
}

static inline void setDirection(Cluster *cluster)
{
    if (directionIsConsistent(cluster))
        cluster->directionFlag = 1;
    else if (directionIsInconsistent(cluster))
        cluster->directionFlag = 2;
}

static inline void setBackbgInfo(Cluster *cluster, Args args)
{
    cluster->bgDiv = args.bgDiv;
    cluster->bgDepth = args.bgDepth;
    cluster->bgReadLen = args.bgReadLen;
}


/*****************
 *** Filtering ***
 *****************/
void intersectBlackList(Cluster *cluster, Args args);


/**********************
 *** Local Assembly ***
 **********************/
#define isLowQualClt(cluster) (((cluster)->isInBlacklist) || ((cluster)->probability <= 0.5))
#define isGermClt(cluster) ((cluster)->cltType == 0)
#define isSomaClt(cluster) ((cluster)->cltType != 0)

int getOuputSegIndex(Cluster *cluster, Segment *segArray, Args args);

#endif // CLUSTER_UTILS_H