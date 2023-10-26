#include "cluster_utils.h"

/**********************
 *** Update Cluster ***
 **********************/
Args initArgs(int numThread, int tid, int minSegLen, int maxDistance, int minOverhang, float bgDiv, float bgDepth, float bgReadLen)
{
    Args args;
    args.numThread = numThread;
    args.tid = tid;
    args.minSegLen = minSegLen;
    args.maxDistance = maxDistance;
    args.minOverhang = minOverhang;
    args.bgDiv = bgDiv;
    args.bgDepth = bgDepth;
    args.bgReadLen = bgReadLen;

    return args;
}

void updateCluster(Cluster *cltArray, Segment *segArray, Args args)
{
    Cluster *cluster = &cltArray[0];

    setTeAlignedFrac(cluster, segArray, args);
    if (cluster->numSeg <= 2) setCltType(cluster, segArray, args);
    if (!isValidCandidate(cluster)) return;

    updateBySegArray(cluster, segArray, args);

    setCltLocationType(cluster, args);
    divideValuesByNumSeg(cluster);
    divideValuesByBackbg(cluster, args);
    setDirection(cluster);
    setBackbgInfo(cluster, args);
}

/************************
 *** Define Candidate ***
 ************************/
void setTeAlignedFrac(Cluster *cluster, Segment *segArray, Args args)
{
    float_t numTeAlignedSeg = 0;

    for (int i = cluster->startIndex; i < cluster->endIndex; i++) {
        if (segArray[i].overhang < args.minOverhang) continue;
        if (segArray[i].numTeAlignment > 0) numTeAlignedSeg += 1;
        cluster->numSeg += 1;
    }

    cluster->teAlignedFrac = numTeAlignedSeg / cluster->numSeg;
}

void setCltType(Cluster *cluster, Segment *segArray, Args args)
{
    if (cluster->numSeg == 1) { cluster->cltType = 1; return; }

    int isFirst = 1;
    for (int i = cluster->startIndex; i < cluster->endIndex; i++) {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang)) continue;
        if (isFirst) {
            bgzf_seek(args.genomeBamFile->fp.bgzf, segment->fileOffset, SEEK_SET);
            bam_read1(args.genomeBamFile->fp.bgzf, args.firstBamRecord);
            isFirst = 0;
            continue;
        }
        bgzf_seek(args.genomeBamFile->fp.bgzf, segment->fileOffset, SEEK_SET);
        bam_read1(args.genomeBamFile->fp.bgzf, args.secondBamRecord);
    }
    
    if (nameIsSame(args.firstBamRecord, args.secondBamRecord)) cluster->cltType = 2;
}

/**********************************
 *** Update Cluster By Segments ***
 **********************************/
void updateBySegArray(Cluster *cluster, Segment *segArray, Args args)
{
    int numLeft = 0, numMiddle = 0, numRight = 0;
    memset(args.teTidCountTable, 0, args.numTeTid * sizeof(int));

    for (int i = cluster->startIndex; i < cluster->endIndex; i++)
    {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang)) continue;
        countValuesFromSeg(cluster, args, segment, &numLeft, &numMiddle, &numRight);
    }

    setEntropy(cluster, numLeft, numMiddle, numRight);
    setBalanceRatio(cluster, numLeft, numRight);
    setNumSegType(cluster);
}

void countValuesFromSeg(Cluster *cluster, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight)
{
    countDifferentSeg(numLeft, numMiddle, numRight, segment);

    cluster->numSegType |= segment->segType;

    if (segment->mapQual < 5) cluster->lowMapQualFrac += 1;
    
    cluster->meanMapQual += segment->mapQual;

    if (isDualClip(segment)) cluster->dualClipFrac += 1;

    countAlnFracs(cluster, segment);

    if (noTeAlignment(segment)) { cluster->meanDivergence += args.bgDiv; return; }

    args.teTidCountTable[segment->teTid] += 1;

    cluster->meanAlnScore += segment->sumAlnScore / segment->numTeAlignment;

    cluster->meanQueryMapFrac += (float_t)segment->sumQueryMapLen / (segment->queryEnd - segment->queryStart);

    cluster->meanDivergence += segment->sumDivergence / segment->numTeAlignment;

    if (directionIsConsistent(segment)) cluster->directionFlag += 1;
    else if (directionIsInconsistent(segment)) cluster->directionFlag += 256;
}

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

static inline void setEntropy(Cluster *cluster, int numLeft, int numMiddle, int numRight)
{ cluster->entropy = getEntropy(numLeft, numMiddle, numRight, cluster->numSeg); }

float_t getEntropy(int numLeft, int numMiddle, int numRight, int numSeg)
{
    float_t entropy = 0;
    if (numLeft > 0)
        entropy -= ((float_t)numLeft / numSeg) * log2((float_t)numLeft / numSeg);
    if (numMiddle > 0)
        entropy -= ((float_t)numMiddle / numSeg) * log2((float_t)numMiddle / numSeg);
    if (numRight > 0)
        entropy -= ((float_t)numRight / numSeg) * log2((float_t)numRight / numSeg);
    return entropy;
}

static inline void setBalanceRatio(Cluster *cluster, int numLeft, int numRight)
{ cluster->balanceRatio = (MIN(numLeft, numRight) + 0.01) / (MAX(numLeft, numRight) + 0.01); }

static inline void setNumSegType(Cluster *cluster)
{ cluster->numSegType = (cluster->numSegType&1) + ((cluster->numSegType&2)>>1) + ((cluster->numSegType&4)>>2); }

/*********************************
 *** Update By Backbg Info ***
 *********************************/
void setCltLocationType(Cluster *cluster, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.repeatAiList, cluster->refStart, cluster->refEnd, 50, &numOverlap, &minDistanceToOverlap);
    ailistQueryInterval(args.gapAiList, cluster->refStart, cluster->refEnd, 50, &numOverlap, &minDistanceToOverlap);

    if (numOverlap == 0) { cluster->locationType = 1; return; }
    if (minDistanceToOverlap < 50) { cluster->locationType = 2; return; }
    cluster->locationType = 4;
}

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
    if (directionIsConsistent(cluster)) cluster->directionFlag = 1;
    else if (directionIsInconsistent(cluster)) cluster->directionFlag = 2;
}

static inline void setBackbgInfo(Cluster *cluster, Args args)
{
    cluster->bgDiv = args.bgDiv;
    cluster->bgDepth = args.bgDepth;
    cluster->bgReadLen = args.bgReadLen;
}