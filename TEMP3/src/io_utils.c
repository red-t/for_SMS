#include "io_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/

/// @brief select a segment to output from cluster
int getOuputSegIdx(Cluster *cluster, Segment *segArray, Args args)
{
    int clipIdx = -1;
    int insIdx = -1;
    int maxReadLen = 0;
    int maxClipLen = 0;

    for (int i = cluster->startIdx; i < cluster->endIdx; i++)
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

    if (insIdx > 0)
        return insIdx;
        
    return clipIdx;
}

/// @brief get extended region of the segment
void setTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize)
{
    *startPtr = 0;
    *endPtr = segment->readLen;

    if (segment->queryStart - flankSize > 0)
        *startPtr = segment->queryStart - flankSize;
    if (segment->queryEnd + flankSize < segment->readLen)
        *endPtr = segment->queryEnd + flankSize;
}

/*************************
 *** Flank Sequence IO ***
 *************************/

FlankRegion initFlankRegion()
{
    FlankRegion region;
    region.start1 = 0;
    region.start2 = 0;
    region.end1 = 0;
    region.end2 = 0;
    return region;
}

/// @brief extract and output flank-seq for all clusters from genome
void extractRefFlanks(char *refFn, Cluster *cltArray, int startIdx, int endIdx)
{
    FlankRegion region = initFlankRegion();
    faidx_t *refFa = fai_load((const char *)refFn);

    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *cluster = &cltArray[i];
        setFlankRegion(cluster, &region);
        outputFlank(cluster, refFa, region);
        outputRefLocal(cluster, refFa, region);
    }

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
}

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *cluster, FlankRegion *region)
{
    // cluster position is a point
    if (cluster->refEnd - cluster->refStart == 1) {
        region->start1 = cluster->refEnd - 450;
        region->end1 = cluster->refEnd + 50;
        region->start2 = cluster->refStart - 50;
        region->end2 = cluster->refStart + 450;
        return;
    }

    // cluster position is a region
    region->start1 = cluster->refEnd - 500;
    region->end1 = cluster->refEnd;
    region->start2 = cluster->refStart;
    region->end2 = cluster->refStart + 500;
    return;
}

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *cluster, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    const char *chrom = faidx_iseq(refFa, cluster->tid);
    char *leftSeq = faidx_fetch_seq64(refFa, chrom, region.start1, region.end1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(refFa, chrom, region.start2, region.end2, &seqLen);

    char *outFn = (char *)malloc(100 * sizeof(char));
    sprintf(outFn, "tmp_anno/%d_%d_flank.fa", cluster->tid, cluster->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);

    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
}

/// @brief output +-500bp around cluster position for tsd annotation
void outputRefLocal(Cluster *cluster, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    int start = (region.start1 < 0) ? 0 : region.start1;
    const char *chrom = faidx_iseq(refFa, cluster->tid);
    char *localSeq = faidx_fetch_seq64(refFa, chrom, start, region.end2, &seqLen);

    char *outFn = (char *)malloc(100 * sizeof(char));
    sprintf(outFn, "tmp_anno/%d_%d_refLocal.fa", cluster->tid, cluster->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">%d\n%s\n", start, localSeq);
    fclose(fp);

    if (localSeq != NULL) {free(localSeq); localSeq=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
}

/*****************************
 *** Insertion Sequence IO ***
 *****************************/

InsRegion initInsRegion()
{
    InsRegion region;
    region.start2 = 0;
    region.end1 = 0;
    region.tid1 = -1;
    region.tid2 = -1;
    region.len1 = 0;
    region.cigar1 = 0;
    region.cigar2 = 0;
    region.flag = 0;
    return region;
}

/// @brief extract and output insertion-seq and flank-seq from assembly
void extractIns(Cluster *cluster)
{
    InsRegion region = initInsRegion();
    setInsRegion(cluster->tid, cluster->idx, &region);
    cluster->flag |= region.flag;
    if (!isFlankMapped(region.flag))
        return;

    char *assmFn = (char *)malloc(100 * sizeof(char));
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", cluster->tid, cluster->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);

    outputInsSeq(cluster, assmFa, region);
    outputTsdSeq(cluster, assmFa, region);

    if (assmFn != NULL) {free(assmFn); assmFn=NULL;}
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
    return;
}

/// @brief define insertion sequence region
void setInsRegion(int cltTid, int cltIdx, InsRegion *region)
{
    char *inputFn = (char *)malloc(100 * sizeof(char));
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
            region->end1 = bam_endpos(bam);
            region->tid1 = bam->core.tid;
            region->len1 = sam_hdr_tid2len(header, region->tid1);
            region->cigar1 = cigarArray[0];
        } else {
            region->start2 = bam->core.pos;
            region->tid2 = bam->core.tid;
            region->cigar2 = cigarArray[numCigar - 1];
        }
    }

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}

    adjustInsRegion(region);
}

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region)
{
    if (region->tid1 < 0 && region->tid2 < 0) // No flank mapped
        return;

    if (isClipInFlank(region->cigar1) && isClipInFlank(region->cigar2)) // Large clip in both flank
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

    if (region->tid1 == region->tid2 && region->end1 <= (region->start2-100)) { // Two flanks map to same contig
        region->flag |= CLT_SAME_FLANK_MAP;
        return;
    }

    if (region->tid1 != region->tid2) // Two flanks map to different contig
        region->flag |= CLT_DIFF_FLANK_MAP;
}

/// @brief output insertion sequence in FASTA format
void outputInsSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region)
{
    char *insSeq = getInsSeq(assmFa, region);
    char *outFn = (char *)malloc(100 * sizeof(char));
    sprintf(outFn, "tmp_anno/%d_%d_insertion.fa", cluster->tid, cluster->idx);

    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0\n%s\n", insSeq);
    fclose(fp);

    if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
}

/// @brief get insertion sequence
char *getInsSeq(faidx_t *assmFa, InsRegion region)
{
    hts_pos_t seqLen;
    if ((region.flag & CLT_LEFT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid1), region.end1, region.len1, &seqLen);

    if ((region.flag & CLT_RIGHT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid2), 0, region.start2, &seqLen);

    if ((region.flag & CLT_SAME_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid1), region.end1, region.start2, &seqLen);

    if ((region.flag & CLT_DIFF_FLANK_MAP) != 0) {
        hts_pos_t len1, len2;
        char *seq1 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid1), region.end1, region.len1, &len1);
        char *seq2 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid2), 0, region.start2, &len2);

        seqLen = len1 + len2 + 100;
        char *temp = (char *)malloc((seqLen + 1) * sizeof(char));
        char *insSeq = temp;
        memcpy(temp, seq1, len1); temp += len1;
        memset(temp, 'N', 100); temp += 100;
        memcpy(temp, seq2, len2);
        temp[seqLen] = '\0';

        if (seq1 != NULL) {free(seq1); seq1=NULL;}
        if (seq2 != NULL) {free(seq2); seq2=NULL;}
        return insSeq;
    }

    return NULL;
}

/// @brief output flank sequences from assembly for tsd annotation
void outputTsdSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region)
{
    if (!isBothFlankMapped(region.flag))
        return;

    hts_pos_t seqLen;
    char *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid1), (region.end1-100), region.end1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, region.tid2), region.start2, (region.start2+100), &seqLen);

    char *outFn = (char *)malloc(100 * sizeof(char));
    sprintf(outFn, "tmp_anno/%d_%d_tsd.fa", cluster->tid, cluster->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);

    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
}