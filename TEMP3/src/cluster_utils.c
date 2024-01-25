#include "cluster_utils.h"

/**********************
 *** Update Cluster ***
 **********************/
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
    float numTeAlignedSeg = 0;

    for (int i = cluster->startIdx; i < cluster->endIdx; i++) {
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
    for (int i = cluster->startIdx; i < cluster->endIdx; i++) {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang)) continue;
        if (isFirst) {
            bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
            bam_read1(args.genomeBam->fp.bgzf, args.firstBamRecord);
            isFirst = 0;
            continue;
        }
        bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
        bam_read1(args.genomeBam->fp.bgzf, args.secondBamRecord);
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

    for (int i = cluster->startIdx; i < cluster->endIdx; i++)
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

    cluster->meanQueryMapFrac += (float)segment->sumQueryMapLen / (segment->queryEnd - segment->queryStart);

    cluster->meanDivergence += segment->sumDivergence / segment->numTeAlignment;

    if (directionIsConsistent(segment)) cluster->directionFlag += 1;
    else if (directionIsInconsistent(segment)) cluster->directionFlag += 256;
}

float getEntropy(int numLeft, int numMiddle, int numRight, int numSeg)
{
    float entropy = 0;
    if (numLeft > 0)
        entropy -= ((float)numLeft / numSeg) * log2((float)numLeft / numSeg);
    if (numMiddle > 0)
        entropy -= ((float)numMiddle / numSeg) * log2((float)numMiddle / numSeg);
    if (numRight > 0)
        entropy -= ((float)numRight / numSeg) * log2((float)numRight / numSeg);
    return entropy;
}


/*********************************
 *** Update By Background Info ***
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


/*****************
 *** Filtering ***
 *****************/
void intersectBlackList(Cluster *cluster, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.blackAiList, cluster->refStart, cluster->refEnd, 2, &numOverlap, &minDistanceToOverlap);

    cluster->isInBlacklist = (uint8_t)(numOverlap > 0);
}
