#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H

#include <stdint.h>
#include <math.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "seg_utils.h"
#include "AIList.h"

/******************************
 *** Cluster related macros ***
 ******************************/
#define CLT_REVERSED            1
#define CLT_IN_BLACKLIST        2
#define CLT_ASSEMBLED           4
#define CLT_LEFT_FLANK_MAP      8
#define CLT_RIGHT_FLANK_MAP     16
#define CLT_DIFF_FLANK_MAP      32
#define CLT_SAME_FLANK_MAP      64
#define CLT_TE_MAP              128
#define CLT_POLYA               256
#define CLT_TSD                 512


/******************
 *** Structures ***
 ******************/

/*! @brief Data container for cluster merged from segments.
 @field  tid                target id of corresponding chromosome
 @field  idx                cluster index in clusters array (0-based)
 @field  refStart           cluster start on reference genome (0-based, included)
 @field  refEnd             cluster end on reference genome (0-based, not-included)
 @field  startIdx           start index in segments array (0-based, include)
 @field  endIdx             end index in segments array (0-based, not-include)
 @field  numSeg             number of segments in the cluster (normalized by bg depth)
 @field  directionFlag      bitwise flag representing cluster direction
                                0: forward
                                1: reverse
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
 @field  flag               bitwise flag representing cluster features
 @field  tsdStart           TSD start on reference genome (0-based, included)
 @field  tsdEnd             TSD end on reference genome (0-based, not-included)
 */
typedef struct Cluster
{
    int         tid;
    int         idx;
    int         refStart;
    int         refEnd;
    int         startIdx;
    int         endIdx;
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
    uint16_t    flag;
    int         tsdStart;
    int         tsdEnd;
} __attribute__((packed)) Cluster;

/// @brief Data container for arguments
typedef struct Args
{
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
    htsFile     *genomeBam;
    bam1_t      *firstBamRecord;
    bam1_t      *secondBamRecord;
    AiList      *repeatAiList;
    AiList      *gapAiList;
    AiList      *blackAiList;
} Args;


/**********************
 *** Update Cluster ***
 **********************/

/// @brief Update cluster values by segments features and background info
void updateCluster(Cluster *cltArray, Segment *segArray, Args args);


/************************
 *** Define Candidate ***
 ************************/
#define overhangIsShort(segment, minOverhang) ((segment)->overhang < (minOverhang))
#define nameIsSame(record1, record2) (strcmp(bam_get_qname((record1)), bam_get_qname((record2))) == 0)
#define isValidCandidate(cluster) ((cluster)->teAlignedFrac >= 0.8)

/// @brief Compute TE-aligned-fraction of a cluster
void setTeAlignedFrac(Cluster *cluster, Segment *segArray, Args args);

/// @brief Set cluster's cltType, which represent the cluster is germ/soma
void setCltType(Cluster *cluster, Segment *segArray, Args args);


/**********************************
 *** Update Cluster By Segments ***
 **********************************/
#define isDualClip(segment) (((segment)->alnType & DUAL_CLIP) == 5)
#define noTeAlignment(segment) ((segment)->numTeAlignment == 0)
#define directionIsConsistent(record) (((record)->directionFlag & 255) > ((record)->directionFlag >> 8))
#define directionIsInconsistent(record) (((record)->directionFlag & 255) < ((record)->directionFlag >> 8))

/// @brief Update cluster values by all segments
void updateBySegArray(Cluster *cluster, Segment *segArray, Args args);

/// @brief Update cluster values by single segment
void countValuesFromSeg(Cluster *cluster, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight);

/// @brief Count number of different type segments
void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment);

/// @brief Count the number of segments with different alnLocationType
void countAlnFracs(Cluster *cluster, Segment *segment);

/// @brief Compute entropy based on the number of different type segments
float getEntropy(int numLeft, int numMiddle, int numRight, int numSeg);

/// @brief Set entropy of the cluster
void setEntropy(Cluster *cluster, int numLeft, int numMiddle, int numRight);

/// @brief Compute BalanceRatio of the cluster
void setBalanceRatio(Cluster *cluster, int numLeft, int numRight);

/// @brief Compute the numer of segment type of the cluster
void setNumSegType(Cluster *cluster);


/*********************************
 *** Update By Background Info ***
 *********************************/

/// @brief Set location type of a cluster
void setCltLocationType(Cluster *cluster, Args args);

/// @brief Divide cluster values by numSeg
void divideByNumSeg(Cluster *cluster);

/// @brief Divide cluster values by background info
void divideByBgInfo(Cluster *cluster, Args args);

/// @brief Set cluster strand
void setDirection(Cluster *cluster);

/// @brief Set background info of a cluster
void setBackbgInfo(Cluster *cluster, Args args);


/*****************
 *** Filtering ***
 *****************/

/// @brief Check if the cluster inersect with blacklist
void intersectBlackList(Cluster *cluster, Args args);


/*******************
 *** Cluster I/O ***
 *******************/
#define isTEMapped(flag) (((flag) & CLT_TE_MAP) != 0)

/// @brief Output formated cluster records
void outputClt(Cluster *cltArray, int startIdx, int endIdx, const char *refFn, const char *teFn);

#endif // CLUSTER_UTILS_H