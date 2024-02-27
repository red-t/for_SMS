#include "cluster_utils.h"

/**********************
 *** Update Cluster ***
 **********************/

/// @brief Update cluster values by segments features and background info
void updateCluster(Cluster *cltArray, Segment *segArray, Args args)
{
    Cluster *cluster = &cltArray[0];
    setTeAlignedFrac(cluster, segArray, args);
    if (cluster->numSeg <= 2)
        setCltType(cluster, segArray, args);
        
    if (!isValidCandidate(cluster))
        return;

    updateBySegArray(cluster, segArray, args);
    setCltLocationType(cluster, args);
    divideByNumSeg(cluster);
    divideByBgInfo(cluster, args);
    setBackbgInfo(cluster, args);
}


/************************
 *** Define Candidate ***
 ************************/

/// @brief Compute TE-aligned-fraction of a cluster
void setTeAlignedFrac(Cluster *cluster, Segment *segArray, Args args)
{
    float numTeAlignedSeg = 0;
    for (int i = cluster->startIdx; i < cluster->endIdx; i++) {
        if (segArray[i].overhang < args.minOverhang)
            continue;
        if (segArray[i].numTeAlignment > 0)
            numTeAlignedSeg += 1;
        cluster->numSeg += 1;
    }

    cluster->teAlignedFrac = numTeAlignedSeg / cluster->numSeg;
}

/// @brief Set cluster's cltType, which represent the cluster is germ/soma
void setCltType(Cluster *cluster, Segment *segArray, Args args)
{
    if (cluster->numSeg == 1) {
        cluster->cltType = 1;
        return;
    }

    int isFirst = 1;
    for (int i = cluster->startIdx; i < cluster->endIdx; i++) {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang))
            continue;
        if (isFirst) {
            bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
            bam_read1(args.genomeBam->fp.bgzf, args.firstBamRecord);
            isFirst = 0;
            continue;
        }
        bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
        bam_read1(args.genomeBam->fp.bgzf, args.secondBamRecord);
    }
    
    if (nameIsSame(args.firstBamRecord, args.secondBamRecord))
        cluster->cltType = 2;
}


/**********************************
 *** Update Cluster By Segments ***
 **********************************/

/// @brief Update cluster values by all segments
void updateBySegArray(Cluster *cluster, Segment *segArray, Args args)
{
    int numLeft = 0, numMiddle = 0, numRight = 0;

    for (int i = cluster->startIdx; i < cluster->endIdx; i++)
    {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang))
            continue;
        countValuesFromSeg(cluster, args, segment, &numLeft, &numMiddle, &numRight);
    }

    setEntropy(cluster, numLeft, numMiddle, numRight);
    setBalanceRatio(cluster, numLeft, numRight);
    setNumSegType(cluster);
}

/// @brief Update cluster values by single segment
void countValuesFromSeg(Cluster *cluster, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight)
{
    countDifferentSeg(numLeft, numMiddle, numRight, segment);
    cluster->numSegType |= segment->segType;
    cluster->meanMapQual += segment->mapQual;
    countAlnFracs(cluster, segment);

    if (segment->mapQual < 5)
        cluster->lowMapQualFrac += 1;
    if (isDualClip(segment))
        cluster->dualClipFrac += 1;
    if (noTeAlignment(segment)) {
        cluster->meanDivergence += args.bgDiv;
        return;
    }

    cluster->meanAlnScore += segment->sumAlnScore / segment->numTeAlignment;
    cluster->meanQueryMapFrac += (float)segment->sumQueryMapLen / (segment->queryEnd - segment->queryStart);
    cluster->meanDivergence += segment->sumDivergence / segment->numTeAlignment;
}

/// @brief Count the number of segments with different type
void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment)
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

/// @brief Count the number of segments with different alnLocationType
void countAlnFracs(Cluster *cluster, Segment *segment)
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

/// @brief Compute entropy based on the number of different type segments
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

/// @brief Set entropy of the cluster
void setEntropy(Cluster *cluster, int numLeft, int numMiddle, int numRight)
{
    cluster->entropy = getEntropy(numLeft, numMiddle, numRight, cluster->numSeg);
    cluster->numSegRaw = cluster->numSeg;
    cluster->numLeft = numLeft;
    cluster->numMiddle = numMiddle;
    cluster->numRight = numRight;
}

/// @brief Compute BalanceRatio of the cluster
void setBalanceRatio(Cluster *cluster, int numLeft, int numRight)
{ cluster->balanceRatio = (MIN(numLeft, numRight) + 0.01) / (MAX(numLeft, numRight) + 0.01); }

/// @brief Compute the numer of segment type of the cluster
void setNumSegType(Cluster *cluster)
{ cluster->numSegType = (cluster->numSegType & 1) + ((cluster->numSegType & 2) >> 1) + ((cluster->numSegType & 4) >> 2); }


/*********************************
 *** Update By Background Info ***
 *********************************/

/// @brief Set location type of a cluster
void setCltLocationType(Cluster *cluster, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.repeatAiList, cluster->refStart, cluster->refEnd, 50, &numOverlap, &minDistanceToOverlap);
    ailistQueryInterval(args.gapAiList, cluster->refStart, cluster->refEnd, 50, &numOverlap, &minDistanceToOverlap);

    if (numOverlap == 0) {
        cluster->locationType = 1;
        return;
    }

    if (minDistanceToOverlap < 50) {
        cluster->locationType = 2;
        return;
    }

    cluster->locationType = 4;
}

/// @brief Divide cluster values by numSeg
void divideByNumSeg(Cluster *cluster)
{
    if (cluster->alnFrac1 > 0)
        cluster->alnFrac1 = cluster->alnFrac1 / cluster->numSeg;
    if (cluster->alnFrac2 > 0)
        cluster->alnFrac2 = cluster->alnFrac2 / cluster->numSeg;
    if (cluster->alnFrac4 > 0)
        cluster->alnFrac4 = cluster->alnFrac4 / cluster->numSeg;
    if (cluster->alnFrac8 > 0)
        cluster->alnFrac8 = cluster->alnFrac8 / cluster->numSeg;
    if (cluster->alnFrac16 > 0)
        cluster->alnFrac16 = cluster->alnFrac16 / cluster->numSeg;

    cluster->dualClipFrac = cluster->dualClipFrac / cluster->numSeg;
    cluster->lowMapQualFrac = cluster->lowMapQualFrac / cluster->numSeg;
    cluster->meanQueryMapFrac = cluster->meanQueryMapFrac / cluster->numSeg;
    cluster->meanDivergence = cluster->meanDivergence / cluster->numSeg;
    cluster->meanMapQual = cluster->meanMapQual / cluster->numSeg;
    cluster->meanAlnScore = cluster->meanAlnScore / cluster->numSeg;
}

/// @brief Divide cluster values by background info
void divideByBgInfo(Cluster *cluster, Args args)
{
    cluster->meanDivergence = cluster->meanDivergence / args.bgDiv;
    cluster->numSeg = cluster->numSeg / args.bgDepth;
}

/// @brief Set background info of a cluster
void setBackbgInfo(Cluster *cluster, Args args)
{
    cluster->bgDiv = args.bgDiv;
    cluster->bgDepth = args.bgDepth;
    cluster->bgReadLen = args.bgReadLen;
}


/*****************
 *** Filtering ***
 *****************/

/// @brief Check if the cluster inersect with blacklist
void intersectBlackList(Cluster *cluster, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.blackAiList, cluster->refStart, cluster->refEnd, 2, &numOverlap, &minDistanceToOverlap);
    if (numOverlap > 0) {
        cluster->isInBlacklist = 1;
        cluster->flag |= CLT_IN_BLACKLIST;
    }
}


/*******************
 *** Cluster I/O ***
 *******************/

/// @brief Output formated cluster records
void outputClt(Cluster *cltArray, int startIdx, int endIdx, const char *refFn, const char *teFn)
{
    faidx_t *refFa = fai_load(refFn);
    faidx_t *teFa = fai_load(teFn);
    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_cltFormated.txt", startIdx);

    char *tsdSeq = NULL, *insSeq = NULL;
    FILE *fp = fopen(outFn, "w");
    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *clt = &cltArray[i];
        if (!isTEMapped(clt->flag))
            continue;

        int isAssembled = ((clt->flag & CLT_ASSEMBLED) != 0);
        tsdSeq = fetchTsdSeq(refFa, clt);
        insSeq = fetchInsSeq(clt);

        fprintf(fp, "%d-%d\t%s\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                clt->tid, clt->idx, faidx_iseq(refFa, clt->tid), clt->refStart, clt->refEnd,
                clt->probability, clt->numSegRaw, clt->numLeft, clt->numMiddle, clt->numRight,
                isAssembled, tsdSeq, insSeq);

        free(tsdSeq);
        free(insSeq);
    }
    fclose(fp);

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
    if (teFa != NULL) {fai_destroy(teFa); teFa=NULL;}
}

/// @brief Fetch TSD sequence from reference genome
char *fetchTsdSeq(faidx_t *refFa, Cluster *clt)
{
    hts_pos_t seqLen;
    char *tsdSeq = NULL;
    if ((clt->flag & CLT_TSD) != 0)
        tsdSeq = faidx_fetch_seq64(refFa, faidx_iseq(refFa, clt->tid), clt->refStart, clt->refEnd-1, &seqLen);
    else {
        tsdSeq = faidx_fetch_seq64(refFa, faidx_iseq(refFa, clt->tid), 0, 0, &seqLen);
        tsdSeq[0] = '.';
    }
    return tsdSeq;
}

/// @brief Fetch insertion sequence from temporary file
char *fetchInsSeq(Cluster *clt)
{
    char insFn[100] = {'\0'};
    sprintf(insFn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);
    faidx_t *insFa = fai_load(insFn);
    
    hts_pos_t seqLen;
    char *insSeq = faidx_fetch_seq64(insFa, faidx_iseq(insFa, 0), 0, INT_MAX, &seqLen);

    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    return insSeq;
}