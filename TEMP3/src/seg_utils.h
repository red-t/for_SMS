#ifndef SEG_UTILS_H
#define SEG_UTILS_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "AIList.h"

/******************
 *** Structures ***
 ******************/
/*! @typedef
 @abstract Structure for segment extracted from CIGAR.
 @field  flag               alignment flag
 @field  mapQual            alignment mapping quality
 @field  queryStart         segment start on query sequence (0-based, include)
 @field  queryEnd           segment end on query sequence (0-based, not-include)
 @field  refPosition        segment position on reference (0-based, include)
 @field  segType            bitwise flag representing segment type
                                1: left clip
                                2: internal insert
                                4: right clip
 @field  alnType            bitwise flag representing read/alignment type, combination of segType(s)
 @field  fileOffset         file offset of the alignment
 @field  alnRefStart        alignment reference start (0-based)
 @field  alnRefEnd          alignment reference end (0-based)
 @field  order              segment order in the same alignment
 @field  numSeg             number of segments from the same alignment
 @field  overhang           length of the shorter anchor part around segment breakpoint
 @field  matchLen           matched length of the alignment (M/=/X)
 @field  readLen            read length of the alignment
 @field  alnLocationType    bitwise flag representing alignment location
                                1:  at least one end inside normal region
                                2:  both ends inside repeat/gap region
                                4:  both ends at repeat/gap boundary
                                8:  one end at repeat/gap boundary, the other inside normal region
                                16: one end at repeat/gap boundary, the other inside repeat/gap
 @field  numTeAlignment     number of alignments mapped to TE
 @field  sumQueryMapLen     summary query mapped-length of the TE alignments (M/I/=/X, no overlap)
 @field  sumAlnScore        summary per base alignment score of the TE alignments
 @field  sumDivergence      summary per base divergence of the TE alignments
 @field  directionFlag      flag representing segment direction
 @field  startIdx           start index in TE alignments array (0-based, include)
 @field  endIdx             end index in TE alignments array (0-based, not-include)
 @field  teTid              majority TE-tid of the segment
 */
typedef struct {
    uint16_t    flag;
    uint8_t     mapQual;
    int         queryStart;
    int         queryEnd;
    int         refPosition;
    uint8_t     segType;
    uint8_t     alnType;
    int64_t     fileOffset;
    int         alnRefStart;
    int         alnRefEnd;
    uint8_t     order;
    uint8_t     numSeg;
    int         overhang;
    int         matchLen;
    int         readLen;
    uint8_t     alnLocationType;
    uint8_t     numTeAlignment;
    int         sumQueryMapLen;
    float       sumAlnScore;
    float       sumDivergence;
    uint16_t    directionFlag;
    int         startIdx;
    int         endIdx;
    int         teTid;
} __attribute__((packed)) Segment;

typedef struct {
    uint8_t     mapQual;
    uint8_t     numSeg;
    uint8_t     segType;
    uint8_t     alnType;
    uint16_t    flag;
    int         queryPosition;
    int         refPosition;
    int         alnRefStart;
    int         matchLen;
    int         readLen;
    int64_t     fileOffset;
} SegValues;

/*! @typedef
 @abstract Structure for TE alignment.
 @field  segIdx         index of corresponding segment in segments array
 @field  AlnScore       alignment score
 @field  queryStart     query start (original direction of segment)
 @field  queryEnd       query end (original direction of segment)
 @field  mapLen         mapped length (M/I/D/=/X)
 @field  divergence     per-base divergence
 @field  flag           bitwise flag of the alignment
 @field  teTid          tid of the TE alignment
 */
typedef struct {
    int     segIdx;
    int     AlnScore;
    int     queryStart;
    int     queryEnd;
    int     mapLen;
    float   divergence;
    int16_t flag;
    int     teTid;
} __attribute__((packed)) TeAlignment;


/***************************
 *** Initialize Segments ***
 ***************************/
#define LEFT_CLIP   1
#define MID_INSERT  2
#define RIGHT_CLIP  4
#define DUAL_CLIP   5

int fillSegmentArray(bam1_t *bamRecord, Segment *segArray, int64_t fileOffset, int minSegLen);
SegValues initSegmentsFromCigar(bam1_t *bamRecord, Segment *segArray, int64_t fileOffset, int minSegLen);
void initSegment(Segment *segment, SegValues segValues, int cigarLen);
void setSameSegValues(Segment *segArray, SegValues segValues);


/***********************
 *** Update Segments ***
 ***********************/
#define isLeftClip(segment) (((segment)->segType & LEFT_CLIP) != 0)
#define isMidInsert(segment) (((segment)->segType & MID_INSERT) != 0)
#define isRightClip(segment) (((segment)->segType & RIGHT_CLIP) != 0)
#define isFirstSegment(segment) ((segment)->order == 0)
#define isSingleSegment(segment) ((segment)->numSeg == 1)
#define isReverse(segment) (((segment)->flag & BAM_FREVERSE) != 0)

void updateSegment(Segment *segArray, AiList *repeatAiList, AiList *gapAiList);
int getOverhang(int overhang, int matchLen);
void setAlnLocationType(AiList *repeatAiList, AiList *gapAiList, Segment *segment);
uint8_t getPointLocationType(AiList *repeatAiList, AiList *gapAiList, int point);
uint8_t getAlnLocationType(uint8_t startLocationType, uint8_t endLocationType);


/***************************************
 *** Update Segments By TeAlignments ***
 ***************************************/
#define isFirstTeAlign(segment) ((segment)->numTeAlignment == 0)
#define isOverlap(teAlignment, prevTeAlignment) ((teAlignment)->queryStart < (prevTeAlignment)->queryEnd)
#define isCover(teAlignment, prevTeAlignment) ((teAlignment)->queryEnd <= (prevTeAlignment)->queryEnd)
#define isSameDirection(segment, teAlignment) (((segment)->flag & BAM_FREVERSE) == ((teAlignment)->flag & BAM_FREVERSE))

void updateSegByTeArray(Segment *segArray, TeAlignment *teArray, int teIdx);

static inline int getQueryMapLen(TeAlignment *teAlignment) { return teAlignment->queryEnd - teAlignment->queryStart; }

static inline int getOverlapQueryMapLen(TeAlignment *teAlignment, TeAlignment *prevTeAlignment)
{
    return teAlignment->queryEnd - prevTeAlignment->queryEnd;
}

void updateSegByTeAlignment(Segment *segment, TeAlignment *teAlignment, int teIdx, int queryMapLen);
void countTeTids(Segment *segArray, TeAlignment *teArray, int *teTidCountTable, int numTeTid);


/*******************************
 *** Initialize TeAlignments ***
 *******************************/
#define alingnedBits 0x3c03f
#define bamIsInvalid(bamRecord) (((bamRecord)->core.flag & BAM_FSECONDARY) != 0 || ((bamRecord)->core.flag & BAM_FUNMAP) != 0)
#define isCigarAligned(cigar) ((alingnedBits >> (bam_cigar_op((cigar))<<1) & 3) & 1)

void fillTeArray(bam1_t *bamRecord, TeAlignment *teArray);
void initQueryPosition(int *queryStartPtr, int *queryEndPtr, bam1_t *bamRecord);

static inline int firstCigarIsClip(uint32_t *cigarArray) { return bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP; }

static inline int lastCigarIsClip(uint32_t *cigarArray, int numCigar) { return bam_cigar_op(cigarArray[numCigar - 1]) == BAM_CSOFT_CLIP; }

void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bamRecord);
float getDivergence(bam1_t *bamRecord, int mapLen);
void initTeAlignment(TeAlignment *teAlignment, bam1_t *bamRecord, int queryStart, int queryEnd, int mapLen, float divergence);


/********************
 *** Trim Segment ***
 ********************/
#define isOdd(number) ((number) & 1)

int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIdx, int sourceStart, int sourceEnd);
int samReallocBamData(bam1_t *bamRecord, size_t desired);

static inline int reallocBamData(bam1_t *bamRecord, size_t desired)
{
    if (desired <= bamRecord->m_data) return 0;
    return samReallocBamData(bamRecord, desired);
}

static inline uint32_t bamGetMemPolicy(bam1_t *bamRecord) { return bamRecord->mempolicy; }

static inline void bamSetMemPolicy(bam1_t *bamRecord, uint32_t policy) { bamRecord->mempolicy = policy; }

void setDestValues(bam1_t *destRecord, int destNameLen, int numNulls, int destDataLen);
uint8_t *setDestName(bam1_t *destRecord, char *destName, int destNameLen, int numNulls);
void copySequence(bam1_t *sourceRecord, bam1_t *destRecord, uint8_t *destDataPtr, int sourceStart, int destSeqLen);


/**********************
 *** Local Assembly ***
 **********************/
void getTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize);

#endif // SEG_UTILS_H