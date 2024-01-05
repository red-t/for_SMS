#include "seg_utils.h"

/***************************
 *** Initialize Segments ***
 ***************************/
int fillSegmentArray(bam1_t *bamRecord, Segment *segArray, int64_t fileOffset, int minSegLen)
{
    SegValues sameSegValues = initSegmentsFromCigar(bamRecord, segArray, fileOffset, minSegLen);
    if (sameSegValues.numSeg > 0) setSameSegValues(segArray, sameSegValues);
    return (int)sameSegValues.numSeg;
}

SegValues initSegmentsFromCigar(bam1_t *bamRecord, Segment *segArray, int64_t fileOffset, int minSegLen)
{
    int numCigar = bamRecord->core.n_cigar;
    int lastCigarIndex = numCigar - 1;
    uint32_t *cigarArray = bam_get_cigar(bamRecord);

    SegValues segValues;
    segValues.numSeg = 0;
    segValues.segType = 0;
    segValues.alnType = 0;
    segValues.matchLen = 0;
    segValues.queryPosition = 0;
    segValues.refPosition = bamRecord->core.pos;

    for (int i = 0; i < numCigar; i++)
    {
        int cigarLen = bam_cigar_oplen(cigarArray[i]);
        switch (bam_cigar_op(cigarArray[i]))
        {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                segValues.queryPosition += cigarLen;
                segValues.refPosition += cigarLen;
                segValues.matchLen += cigarLen;
                break;

            case BAM_CSOFT_CLIP:
            case BAM_CINS:
                if (cigarLen < minSegLen) { segValues.queryPosition += cigarLen; break; }

                if (i == 0) segValues.segType = LEFT_CLIP;
                else if (i == lastCigarIndex) segValues.segType = RIGHT_CLIP;
                else segValues.segType = MID_INSERT;
                
                initSegment(&segArray[segValues.numSeg], segValues, cigarLen);
                segValues.numSeg += 1;
                segValues.alnType |= segValues.segType;
                segValues.queryPosition += cigarLen;
                break;

            case BAM_CDEL:
            case BAM_CREF_SKIP:
                segValues.refPosition += cigarLen;
                break;
            default:
                break;
        }
    }

    segValues.fileOffset = fileOffset;
    segValues.mapQual = bamRecord->core.qual;
    segValues.flag = bamRecord->core.flag;
    segValues.alnRefStart = bamRecord->core.pos;
    segValues.readLen = bamRecord->core.l_qseq;
    return segValues;
}

void initSegment(Segment *segment, SegValues segValues, int cigarLen)
{
    segment->queryStart = segValues.queryPosition;
    segment->queryEnd = segValues.queryPosition + cigarLen;
    segment->refPosition = segValues.refPosition;
    segment->overhang = segValues.matchLen;
    segment->segType = segValues.segType;
    segment->order = segValues.numSeg;
}

void setSameSegValues(Segment *segArray, SegValues segValues)
{
    for (uint8_t i = 0; i < segValues.numSeg; i++)
    {
        Segment *segment = &segArray[i];
        segment->mapQual = segValues.mapQual;
        segment->flag = segValues.flag;
        segment->alnType = segValues.alnType;
        segment->numSeg = segValues.numSeg;
        segment->alnRefStart = segValues.alnRefStart;
        segment->alnRefEnd = segValues.refPosition;
        segment->matchLen = segValues.matchLen;
        segment->readLen = segValues.readLen;
        segment->fileOffset = segValues.fileOffset;
    }
}


/***********************
 *** Update Segments ***
 ***********************/
void updateSegment(Segment *segArray, AiList *repeatAiList, AiList *gapAiList)
{
    Segment *segment = &segArray[0];

    if (isLeftClip(segment))
        segment->overhang = segment->matchLen; // update overhang
    else if (isMidInsert(segment))
        segment->overhang = getOverhang(segment->overhang, segment->matchLen);

    if (!isFirstSegment(segment)) return;
    setAlnLocationType(repeatAiList, gapAiList, segment); // update alnLocationType
    
    if (isSingleSegment(segment)) return;
    for (uint8_t i = 1; i < segment->numSeg; i++) // update alnLocationType for other segments
        segArray[i].alnLocationType = segment->alnLocationType;
}

int getOverhang(int overhang, int matchLen)
{ return (matchLen - overhang) < overhang ? (matchLen - overhang) : overhang; }

void setAlnLocationType(AiList *repeatAiList, AiList *gapAiList, Segment *segment)
{
    uint8_t startLocationType = getPointLocationType(repeatAiList, gapAiList, segment->alnRefStart);
    uint8_t endLocationType = getPointLocationType(repeatAiList, gapAiList, segment->alnRefEnd);
    segment->alnLocationType = getAlnLocationType(startLocationType, endLocationType);
}

uint8_t getPointLocationType(AiList *repeatAiList, AiList *gapAiList, int point)
{
    int numOverlap = 0;
    int minDistanceToOverlap = 0x7fffffff;
    ailistQueryPoint(repeatAiList, point, 50, &numOverlap, &minDistanceToOverlap);
    ailistQueryPoint(gapAiList, point, 50, &numOverlap, &minDistanceToOverlap);

    if (numOverlap == 0) return 1;
    if (minDistanceToOverlap < 50) return 2;
    return 4;
}

uint8_t getAlnLocationType(uint8_t startLocationType, uint8_t endLocationType)
{
    switch (startLocationType | endLocationType)
    {
        case 1:
        case 5:
            return 1;   // at least one end inside normal region
        case 4:
            return 2;   // both ends inside repeat/gap region
        case 2:
            return 4;   // both ends at repeat/gap boundary
        case 3:
            return 8;   // one end at repeat/gap boundary, the other inside normal region
        case 6:
            return 16;  // one end at repeat/gap boundary, the other inside repeat/gap
        default:
            return 0;
    }
}


/***************************************
 *** Update Segments By TeAlignments ***
 ***************************************/
void updateSegByTeArray(Segment *segArray, TeAlignment *teArray, int teIndex)
{
    TeAlignment *teAlignment = &teArray[teIndex];
    int segIndex = teAlignment->segIndex;
    Segment *segment = &segArray[segIndex];

    if (isFirstTeAlign(segment)) {
        segment->startIndex = teIndex;
        int queryMapLen = getQueryMapLen(teAlignment);
        updateSegByTeAlignment(segment, teAlignment, teIndex, queryMapLen);
        return;
    }

    int prevIndex = teIndex - 1;
    TeAlignment *prevTeAlignment = &teArray[prevIndex];
    
    if (!isOverlap(teAlignment, prevTeAlignment)) {
        int queryMapLen = getQueryMapLen(teAlignment);
        updateSegByTeAlignment(segment, teAlignment, teIndex, queryMapLen);
        return;
    }

    if (isCover(teAlignment, prevTeAlignment)) {
        teAlignment->queryEnd = prevTeAlignment->queryEnd;
        return;
    }

    int queryMapLen = getOverlapQueryMapLen(teAlignment, prevTeAlignment);
    updateSegByTeAlignment(segment, teAlignment, teIndex, queryMapLen);
}

static inline int getQueryMapLen(TeAlignment *teAlignment)
{ return teAlignment->queryEnd - teAlignment->queryStart; }

static inline int getOverlapQueryMapLen(TeAlignment *teAlignment, TeAlignment *prevTeAlignment)
{ return teAlignment->queryEnd - prevTeAlignment->queryEnd; }

void updateSegByTeAlignment(Segment *segment, TeAlignment *teAlignment, int teIndex, int queryMapLen)
{
    segment->endIndex = teIndex + 1;
    segment->numTeAlignment += 1;
    segment->sumQueryMapLen += queryMapLen;
    segment->sumDivergence += teAlignment->divergence;
    segment->sumAlnScore += (float)teAlignment->AlnScore / teAlignment->mapLen;

    if (isReverse(segment))
    {
        if (isSameDirection(segment, teAlignment)) { segment->directionFlag += (1 << 8); return; }
        segment->directionFlag += 1; return;
    }
    if (isSameDirection(segment, teAlignment)) { segment->directionFlag += 1; return; }
    segment->directionFlag += (1 << 8);
}

void countTeTids(Segment *segArray, TeAlignment *teArray, int *teTidCountTable, int numTeTid)
{
    Segment *segment = &segArray[0];
    memset(teTidCountTable, 0, numTeTid * sizeof(int));

    for (int i = segment->startIndex; i < segment->endIndex; i++)
        teTidCountTable[teArray[i].teTid] += 1;
}


/*******************************
 *** Initialize TeAlignments ***
 *******************************/
void fillTeArray(bam1_t *bamRecord, TeAlignment *teArray)
{
    int queryStart, queryEnd;
    initQueryPosition(&queryStart, &queryEnd, bamRecord);

    int mapLen;
    float divergence;
    getMapLenAndDiv(&mapLen, &divergence, bamRecord);

    initTeAlignment(&teArray[0], bamRecord, queryStart, queryEnd, mapLen, divergence);
}

void initQueryPosition(int *queryStartPtr, int *queryEndPtr, bam1_t *bamRecord)
{
    int numCigar = bamRecord->core.n_cigar;
    uint32_t *cigarArray = bam_get_cigar(bamRecord);
    *queryStartPtr = 0;
    *queryEndPtr = bamRecord->core.l_qseq;

    if (!bam_is_rev(bamRecord)) {
        if (firstCigarIsClip(cigarArray))
            *queryStartPtr = bam_cigar_oplen(cigarArray[0]);
        if (lastCigarIsClip(cigarArray, numCigar))
            *queryEndPtr = bamRecord->core.l_qseq - bam_cigar_oplen(cigarArray[numCigar-1]);
        return;
    }

    if (lastCigarIsClip(cigarArray, numCigar))
        *queryStartPtr = bam_cigar_oplen(cigarArray[numCigar-1]);
    if (firstCigarIsClip(cigarArray))
        *queryEndPtr = bamRecord->core.l_qseq - bam_cigar_oplen(cigarArray[0]);
}

static inline int firstCigarIsClip(uint32_t *cigarArray)
{ return bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP; }

static inline int lastCigarIsClip(uint32_t *cigarArray, int numCigar)
{ return bam_cigar_op(cigarArray[numCigar-1]) == BAM_CSOFT_CLIP; }

void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bamRecord)
{
    uint32_t *cigarArray = bam_get_cigar(bamRecord);
    int numCigar = bamRecord->core.n_cigar;
    int i, mapLen;

    for (i=mapLen=0; i < numCigar; i++)
        if (isCigarAligned(cigarArray[i])) mapLen += bam_cigar_oplen(cigarArray[i]);

    *mapLenPtr = mapLen;
    *divergencePtr = getDivergence(bamRecord, mapLen);
}

float getDivergence(bam1_t *bamRecord, int mapLen)
{ return (float)(bam_aux2i(bam_aux_get(bamRecord, "NM")) - bam_aux2i(bam_aux_get(bamRecord, "nn"))) / mapLen; }

void initTeAlignment(TeAlignment *teAlignment, bam1_t *bamRecord, int queryStart, int queryEnd, int mapLen, float divergence)
{
    teAlignment->segIndex = atoi(bam_get_qname(bamRecord));
    teAlignment->AlnScore = bam_aux2i(bam_aux_get(bamRecord, "AS"));
    teAlignment->queryStart = queryStart;
    teAlignment->queryEnd = queryEnd;
    teAlignment->mapLen = mapLen;
    teAlignment->divergence = divergence;
    teAlignment->flag = bamRecord->core.flag;
    teAlignment->teTid = bamRecord->core.tid;
}


/********************
 *** Trim Segment ***
 ********************/
int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIndex, int sourceStart, int sourceEnd)
{
    // use segIndex as destName
    char destName[12];
    sprintf(destName, "%d", segIndex);
    
    int destNameLen = strlen(destName);
    int numNulls = 4 - destNameLen % 4;
    int destSeqLen = sourceEnd - sourceStart;
    int destDataLen = destNameLen + numNulls + (destSeqLen + 1)/2 + 1;

    // re-allocate the data buffer as needed.
    if (reallocBamData(destRecord, destDataLen) < 0) return -1;

    setDestValues(destRecord, destNameLen, numNulls, destDataLen);
    uint8_t *destDataPtr = setDestName(destRecord, destName, destNameLen, numNulls);
    // if sourceRecord is reverse, the sequence is reverse-complementary
    copySequence(sourceRecord, destRecord, destDataPtr, sourceStart, destSeqLen);

    return (int)destDataLen;
}

static inline int reallocBamData(bam1_t *bamRecord, size_t desired)
{
    if (desired <= bamRecord->m_data) return 0;
    return samReallocBamData(bamRecord, desired);
}

int samReallocBamData(bam1_t *bamRecord, size_t desired)
{
    // similar to sam_realloc_bam_data in htslib sam.c
    uint32_t newDataMem;
    uint8_t *newData;
    newDataMem = desired;
    kroundup32(newDataMem);
    if (newDataMem < desired) { errno = ENOMEM; return -1; }
    if ((bamGetMemPolicy(bamRecord) & BAM_USER_OWNS_DATA) == 0) {
        newData = realloc(bamRecord->data, newDataMem);
    } else {
        if ((newData = malloc(newDataMem)) != NULL) {
            if (bamRecord->l_data > 0)
                memcpy(newData, bamRecord->data,
                       bamRecord->l_data < (int)bamRecord->m_data ? bamRecord->l_data : (int)bamRecord->m_data);
            bamSetMemPolicy(bamRecord, bamGetMemPolicy(bamRecord) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!newData) return -1;
    bamRecord->data = newData;
    bamRecord->m_data = newDataMem;
    return 0;
}

static inline uint32_t bamGetMemPolicy(bam1_t *bamRecord)
{ return bamRecord->mempolicy; }

static inline void bamSetMemPolicy(bam1_t *bamRecord, uint32_t policy)
{ bamRecord->mempolicy = policy; }

void setDestValues(bam1_t *destRecord, int destNameLen, int numNulls, int destDataLen)
{
    destRecord->core.l_qname = (uint16_t)(destNameLen + numNulls);
    destRecord->core.l_extranul = (uint8_t)(numNulls - 1);
    destRecord->l_data = (int)destDataLen;
    destRecord->core.flag = BAM_FUNMAP;
}

uint8_t *setDestName(bam1_t *destRecord, char *destName, int destNameLen, int numNulls)
{
    uint8_t *destDataPtr = destRecord->data;
    
    memcpy(destDataPtr, destName, destNameLen);
    for (int i = 0; i < numNulls; i++)
        destDataPtr[destNameLen + i] = '\0';
    
    destDataPtr += destNameLen + numNulls;
    return destDataPtr;
}

void copySequence(bam1_t *sourceRecord, bam1_t *destRecord, uint8_t *destDataPtr, int sourceStart, int destSeqLen)
{
    uint8_t *sourceSeqPtr = bam_get_seq(sourceRecord);
    sourceSeqPtr += (sourceStart >> 1);

    // 1. when sourceStart is even, segment actually start from 1-st base of sourceSeqPtr[0].
    //    else, start from 2-nd base of sourceSeqPtr[0].
    // 2. when destSeqLen is even, memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1) will copy destSeqLen bases.
    //    else, will copy (destSeqLen+1) bases.
    if (!isOdd(sourceStart)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
        destRecord->core.l_qseq = (int)destSeqLen;
        return;
    }

    if (!isOdd(destSeqLen)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+2) >> 1);
        destRecord->core.l_qseq = (int)destSeqLen + 1;
        return;
    }

    memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
    destRecord->core.l_qseq = (int)destSeqLen + 1;
}

/**********************
 *** Local Assembly ***
 **********************/
void getTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize)
{
    *startPtr = 0;
    *endPtr = segment->readLen;

    if (segment->queryStart - flankSize > 0)
        *startPtr = segment->queryStart - flankSize;
    if (segment->queryEnd + flankSize < segment->readLen)
        *endPtr = segment->queryEnd + flankSize;
}