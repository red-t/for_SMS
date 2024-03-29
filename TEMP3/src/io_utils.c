#include "io_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/

/// @brief select a segment to output as assembled contig
int getOuputSegIdx(Cluster *clt, Segment *segArray, Args args)
{
    int clipIdx = -1;
    int insIdx = -1;
    int maxReadLen = 0;
    int maxClipLen = 0;

    for (int i = clt->startIdx; i < clt->endIdx; i++)
    {
        Segment *segment = &segArray[i];
        if (overhangIsShort(segment, args.minOverhang))
            continue;

        if (isMidInsert(segment)) {
            if (segment->readLen <= maxReadLen)
                continue;

            maxReadLen = segment->readLen;
            insIdx = i;
            continue;
        }

        int clipLen = segment->queryEnd - segment->queryStart;
        if (clipLen > maxClipLen) {
            maxClipLen = clipLen;
            clipIdx = i;
        }
    }

    if (insIdx > 0) // select longest spanning read
        return insIdx;

    return clipIdx; // or clip read with largest clipSize
}

/// @brief get extended region of the segment
void setTrimRegion(Segment *segment, int *start, int *end, int flankSize)
{
    *start = 0;
    *end = segment->readLen;

    if (segment->queryStart - flankSize > 0)
        *start = segment->queryStart - flankSize;
    if (segment->queryEnd + flankSize < segment->readLen)
        *end = segment->queryEnd + flankSize;
}


/*************************
 *** Flank Sequence IO ***
 *************************/

/// @brief Initiate FlankRegion
FlankRegion initFlankRegion()
{
    FlankRegion region;
    region.start1 = 0;
    region.start2 = 0;
    region.end1 = 0;
    region.end2 = 0;
    return region;
}

/// @brief output flank-seq and local-seq for all clusters
void extractRefFlanks(char *refFn, Cluster *cltArray, int startIdx, int endIdx)
{
    FlankRegion region = initFlankRegion();
    faidx_t *refFa = fai_load((const char *)refFn);

    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *clt = &cltArray[i];
        setFlankRegion(clt, &region);
        outputFlank(clt, refFa, region);
        outputLocal(clt, refFa, region);
    }

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
}

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *clt, FlankRegion *region)
{
    region->start1 = clt->refStart - 499;
    region->end1 = clt->refStart;
    region->start2 = clt->refEnd;
    region->end2 = clt->refEnd + 499;
    return;
}

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *clt, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    const char *chrom = faidx_iseq(refFa, clt->tid);
    char *leftSeq = faidx_fetch_seq64(refFa, chrom, region.start1, region.end1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(refFa, chrom, region.start2, region.end2, &seqLen);

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_flank.fa", clt->tid, clt->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);

    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
}

/// @brief output +-500bp local-seq around cluster position for tsd annotation
void outputLocal(Cluster *clt, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    int start = (region.start1 < 0) ? 0 : region.start1;
    const char *chrom = faidx_iseq(refFa, clt->tid);
    char *localSeq = faidx_fetch_seq64(refFa, chrom, start, region.end2, &seqLen);

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_local.fa", clt->tid, clt->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">%d\n%s\n", start, localSeq);
    fclose(fp);

    if (localSeq != NULL) {free(localSeq); localSeq=NULL;}
}


/*****************************
 *** Insertion Sequence IO ***
 *****************************/

/// @brief Initiate InsRegion
InsRegion initInsRegion()
{
    InsRegion region;
    region.rightMost = 0;
    region.leftMost = 0;
    region.tid1 = -1;
    region.tid2 = -1;
    region.cigar1 = 0;
    region.cigar2 = 0;
    region.flag = 0;
    return region;
}

/// @brief output insertion-seq and tsd-containing-seq from contig
void extractIns(Cluster *clt)
{
    InsRegion region = initInsRegion();
    setInsRegion(clt->tid, clt->idx, &region);
    clt->flag |= region.flag;
    if (!isFlankMapped(region.flag))
        return;

    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);
    outputInsSeq(clt, assmFa, region);

    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
    return;
}

/// @brief define insertion-seq region by Flank-To-Assm alignments
void setInsRegion(int cltTid, int cltIdx, InsRegion *region)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_FlankToAssm.bam", cltTid, cltIdx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam) || bamIsSup(bam) || bam_is_rev(bam))
            continue;

        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArray = bam_get_cigar(bam);

        if (isLeftFlank(bam)) {
            region->leftMost = bam_endpos(bam);
            region->tid1 = bam->core.tid;
            region->cigar1 = cigarArray[0];
        } else {
            region->rightMost = bam->core.pos;
            region->tid2 = bam->core.tid;
            region->cigar2 = cigarArray[numCigar - 1];
        }
    }

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}

    adjustInsRegion(region);
}

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region)
{
    if (region->tid1 < 0 && region->tid2 < 0) // No flank mapped
        return;

    // Large clip in both flank
    if (isClipInFlank(region->cigar1) && isClipInFlank(region->cigar2))
        return;

    if (region->tid2 < 0) { // Only left lank mapped
        region->flag |= CLT_LEFT_FLANK_MAP;
        return;
    }

    if (region->tid1 < 0) { // Only right lank mapped
        region->flag |= CLT_RIGHT_FLANK_MAP;
        return;
    }

    if (isClipInFlank(region->cigar2)) { // Large flank in right flank
        region->flag |= CLT_LEFT_FLANK_MAP;
        return;
    }

    if (isClipInFlank(region->cigar1)) { // Large flank in left flank
        region->flag |= CLT_RIGHT_FLANK_MAP;
        return;
    }

    // Two flanks map to same contig
    if (region->tid1 == region->tid2 && region->leftMost <= (region->rightMost-100)) {
        region->flag |= CLT_SAME_FLANK_MAP;
        return;
    }

    if (region->tid1 != region->tid2) // Two flanks map to different contig
        region->flag |= CLT_DIFF_FLANK_MAP;
}

/// @brief output insertion-seq for annotation
void outputInsSeq(Cluster *clt, faidx_t *assmFa, InsRegion region)
{
    clt->tid1 = region.tid1; clt->leftMost = region.leftMost;
    clt->tid2 = region.tid2; clt->rightMost = region.rightMost;
    char *insSeq = getInsSeq(assmFa, clt);
    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);

    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">%d_%d_%d_%d\n%s\n", region.tid1, region.leftMost, region.tid2, region.rightMost, insSeq);
    fclose(fp);

    if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
}

/// @brief get insertion-seq
char *getInsSeq(faidx_t *assmFa, Cluster *clt)
{
    hts_pos_t seqLen;
    if ((clt->flag & CLT_LEFT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, INT_MAX, &seqLen);

    if ((clt->flag & CLT_RIGHT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, (clt->rightMost-1), &seqLen);

    if ((clt->flag & CLT_SAME_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, (clt->rightMost-1), &seqLen);

    if ((clt->flag & CLT_DIFF_FLANK_MAP) != 0) {
        char *seq1 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, INT_MAX, &seqLen);
        char *seq2 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, (clt->rightMost-1), &seqLen);

        int len1 = strlen(seq1), len2 = strlen(seq2);
        seqLen = len1 + len2 + 500;
        char *insSeq = (char *)malloc((seqLen + 1) * sizeof(char));
        char *temp = insSeq;
        memcpy(temp, seq1, len1); temp += len1;
        memset(temp, 'N', 500); temp += 500;
        memcpy(temp, seq2, len2);
        insSeq[seqLen] = '\0';

        if (seq1 != NULL) {free(seq1); seq1=NULL;}
        if (seq2 != NULL) {free(seq2); seq2=NULL;}
        return insSeq;
    }

    return NULL;
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

    int isAssembled;
    char *tsdSeq = NULL, *insSeq = NULL, *leftSeq = NULL, *rightSeq = NULL;
    FILE *fp = fopen(outFn, "w");
    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *clt = &cltArray[i];
        if (!isTEMapped(clt->flag))
            continue;

        isAssembled = ((clt->flag & CLT_ASSEMBLED) != 0);
        tsdSeq = fetchTsdSeq(refFa, clt);
        fetchSeqs(clt, &insSeq, &leftSeq, &rightSeq);

        fprintf(fp, "%d-%d\t%s\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\n",
                clt->tid, clt->idx, faidx_iseq(refFa, clt->tid), clt->refStart, clt->refEnd,
                clt->probability, clt->numSegRaw, clt->numLeft, clt->numMiddle, clt->numRight,
                isAssembled, tsdSeq, insSeq, leftSeq, rightSeq, clt->flag);

        free(tsdSeq); free(insSeq); free(leftSeq); free(rightSeq);
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

/// @brief Fetch flank sequence from temporary file
void fetchSeqs(Cluster *clt, char **insSeq, char **leftSeq, char **rightSeq)
{
    hts_pos_t seqLen;
    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);

    if ((clt->flag & CLT_LEFT_FLANK_MAP) != 0) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        **rightSeq = '.';
        goto END;
    }

    if ((clt->flag & CLT_RIGHT_FLANK_MAP) != 0) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
        **leftSeq = '.';
        goto END;
    }
    
    *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
    *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
    goto END;

    END:
    *insSeq = getInsSeq(assmFa, clt);
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}