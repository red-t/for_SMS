#include "io_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/

/// @brief select a segment to output as assembled contig
int getOuputSegIdx(Cluster *clt, Segment *segArr, Args args)
{
    int clipIdx = -1;
    int insIdx = -1;
    int maxReadLen = 0;
    int maxClipLen = 0;

    for (int i = clt->startIdx; i < clt->endIdx; i++)
    {
        Segment *segment = &segArr[i];
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
void extractRefFlanks(char *refFn, Cluster *cltArr, int startIdx, int endIdx)
{
    FlankRegion region = initFlankRegion();
    faidx_t *refFa = fai_load((const char *)refFn);

    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *clt = &cltArr[i];
        setFlankRegion(clt, &region);
        outputFlank(clt, refFa, region);
        // outputLocal(clt, refFa, region);
        outputLocal(clt, refFa, 200, "local");

        // Output local region for secodary defining
        outputLocal(clt, refFa, 2000, "ref");
        
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

// /// @brief output +-500bp local-seq around cluster position for tsd annotation
// void outputLocal(Cluster *clt, faidx_t *refFa, FlankRegion region)
// {
//     hts_pos_t seqLen;
//     int start = (region.start1 < 0) ? 0 : region.start1;
//     const char *chrom = faidx_iseq(refFa, clt->tid);
//     char *localSeq = faidx_fetch_seq64(refFa, chrom, start, region.end2, &seqLen);

//     char outFn[100] = {'\0'};
//     sprintf(outFn, "tmp_anno/%d_%d_local.fa", clt->tid, clt->idx);
//     FILE *fp = fopen(outFn, "w");
//     fprintf(fp, ">%d\n%s\n", start, localSeq);
//     fclose(fp);

//     if (localSeq != NULL) {free(localSeq); localSeq=NULL;}
// }

/// @brief Output local region around refrence breakpoint with length of 2*length
void outputLocal(Cluster *clt, faidx_t *refFa, int length, const char *suffix)
{
    hts_pos_t seqLen;
    int start = (clt->refStart - length < 0) ? 0 : clt->refStart - length;
    int end = clt->refStart + length - 1;
    const char *chrom = faidx_iseq(refFa, clt->tid);
    char *localSeq = faidx_fetch_seq64(refFa, chrom, start, end, &seqLen);

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_%s.fa", clt->tid, clt->idx, suffix);
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

// /// @brief output insertion-seq from contig
// void extractIns(Cluster *clt)
// {
//     InsRegion region = initInsRegion();
//     setInsRegion(clt, &region);
//     if (!isFlankMapped(clt->flag))
//         return;

//     char assmFn[100] = {'\0'};
//     sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
//     faidx_t *assmFa = fai_load((const char *)assmFn);

//     if (isBothFlankMapped(clt->flag))
//         outputInsSeq(assmFa, clt);
//     else
//         outputAssmFlank(assmFa, clt);

//     if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
//     return;
// }

/// @brief Primary defining of insertion region
void defineInsRegion(char *refFn, Cluster *clt)
{
    // Primary defining
    InsRegion region = initInsRegion();
    setInsRegion(clt, &region);
}

/// @brief define insertion-seq region by Flank-To-Assm alignments
void setInsRegion(Cluster *clt, InsRegion *region)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_FlankToAssm.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int minRightClip = INT_MAX;
    int minLeftClip = INT_MAX;
    int clipSize = 0;
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam) || bam_is_rev(bam))
            continue;

        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArr = bam_get_cigar(bam);

        if (isLeftFlank(bam)) {
            clipSize = isClipInFlank(cigarArr[numCigar - 1], 0) ? bam_cigar_oplen(cigarArr[numCigar - 1]) : 0;
            if (clipSize > minRightClip)
                continue;
            minRightClip = clipSize;
            region->leftMost = bam_endpos(bam);
            region->tid1 = bam->core.tid;
            region->cigar1 = cigarArr[0];
        } else {
            clipSize = isClipInFlank(cigarArr[0], 0) ? bam_cigar_oplen(cigarArr[0]) : 0;
            if (clipSize > minLeftClip)
                continue;
            minLeftClip = clipSize;
            bam_cigar_oplen(cigarArr[0]);
            region->rightMost = bam->core.pos;
            region->tid2 = bam->core.tid;
            region->cigar2 = cigarArr[numCigar - 1];
        }
    }

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}

    adjustInsRegion(clt, region);
}

/// @brief Adjust insertion region and clt->flag
void adjustInsRegion(Cluster *clt, InsRegion *region)
{
    clt->tid1 = region->tid1;
    clt->leftMost = region->leftMost;
    clt->tid2 = region->tid2;
    clt->rightMost = region->rightMost;

    // Only left end can be defined
    if ((region->tid2 < 0 || isClipInFlank(region->cigar2, 400)) && (region->tid1 >= 0 && !isClipInFlank(region->cigar1, 400))) {
        clt->flag |= CLT_LEFT_FLANK_MAP;
        return;
    }

    // Only right end can be defined
    if ((region->tid1 < 0 || isClipInFlank(region->cigar1, 400)) && (region->tid2 >= 0 && !isClipInFlank(region->cigar2, 400))) {
        clt->flag |= CLT_RIGHT_FLANK_MAP;
        return;
    }

    // Both ends can be defined
    if ((region->tid1 >= 0 && region->tid2 >= 0) && (!isClipInFlank(region->cigar1, 400) && !isClipInFlank(region->cigar2, 400))) {
        if (region->leftMost > region->rightMost - 100)
            return;
        clt->flag |= (region->tid1 == region->tid2) ? CLT_SAME_FLANK_MAP : CLT_DIFF_FLANK_MAP;
        return;
    }
}

/// @brief Secondary defining of insertion region
void refineInsRegion(Cluster *clt)
{
    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);
    
    reSetInsRegion1(clt, assmFa);
    if (isFlankMapped(clt->flag))
        outputInsSeq(assmFa, clt);

    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}

/// @brief Reset insertion region
void reSetInsRegion1(Cluster *clt, faidx_t *assmFa)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_AssmToRef.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    // Traverse CIGAR to detect insertion site for each aligment
    int localStart = atoi(sam_hdr_tid2name(header, 0));
    int refLen = header->target_len[0];
    if (refLen == 4000) {
        clt->refStart = localStart + 2000;
        clt->refEnd = clt->refStart + 1;
    }

    int tid1 = -1, refPos1 = 0, leftMost = 0;
    int tid2 = -1, refPos2 = 0, rightMost = 0;

    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsUnmap(bam) || bam_is_rev(bam))
            continue;

        int queryPosition = 0;                        // position on assembly
        int refPosition = localStart + bam->core.pos; // position on local reference
        int numCigar = bam->core.n_cigar;
        int lastCigarIdx = numCigar - 1;
        uint32_t *cigarArr = bam_get_cigar(bam);
        for (int i = 0; i < numCigar; i++)
        {
            int cigarLen = bam_cigar_oplen(cigarArr[i]);
            switch (bam_cigar_op(cigarArr[i]))
            {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                queryPosition += cigarLen;
                refPosition += cigarLen;
                break;

            case BAM_CSOFT_CLIP:
            case BAM_CINS:
                if (cigarLen < 100)
                    goto SKIP;

                int deltaDistToBp1 = abs(refPosition - clt->refStart) - abs(refPos1 - clt->refStart);
                int deltaDistToBp2 = abs(refPosition - clt->refStart) - abs(refPos2 - clt->refStart);

                if (i == 0) {
                    // left-clip, update rightMost
                    if (deltaDistToBp2 > 0)
                        goto SKIP;

                    rightMost = cigarLen;
                    refPos2 = refPosition;
                    tid2 = getTid(bam, assmFa);
                } else if (i == lastCigarIdx) {
                    // right-clip, update leftMost
                    if (deltaDistToBp1 > 0)
                        goto SKIP;

                    leftMost = queryPosition;
                    refPos1 = refPosition;
                    tid1 = getTid(bam, assmFa);
                } else {
                    // insertion, update leftMost && rightMost
                    if (deltaDistToBp1 > 0 || deltaDistToBp2 > 0)
                        goto SKIP;
                    
                    leftMost = queryPosition;
                    refPos1 = refPosition;
                    tid1 = getTid(bam, assmFa);

                    rightMost = queryPosition + cigarLen;
                    refPos2 = refPosition;
                    tid2 = getTid(bam, assmFa);
                }

                SKIP:
                    queryPosition += cigarLen;
                    break;

            case BAM_CDEL:
            case BAM_CREF_SKIP:
                refPosition += cigarLen;
                break;
                
            default:
                break;
            }
        }
    }
    
    reAdjustInsRegion(clt, tid1, refPos1, leftMost, tid2, refPos2, rightMost);

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
}

/// @brief Get tid in the assmFa, based on qname of the bam record
int getTid(bam1_t *bam, faidx_t *assmFa)
{
    int nseq = faidx_nseq(assmFa);
    const char *seqName = bam_get_qname(bam);
    int tid = -1;

    for (int i = 0; i < nseq; i++)
    {
        const char *currentName = faidx_iseq(assmFa, i);
        if (strcmp(currentName, seqName) == 0) {
            tid = i;
            break;
        }
    }

    return tid;
}

/// @brief Re-adjust insertion region and clt->flag
void reAdjustInsRegion(Cluster *clt, int tid1, int refPos1, int leftMost, int tid2, int refPos2, int rightMost)
{
    int distToBp1 = abs(refPos1 - clt->refStart);
    int distToBp2 = abs(refPos2 - clt->refStart);

    // Only left end can be refined
    if ((tid2 < 0 || distToBp2 > 50) && (tid1 >= 0 && distToBp1 <= 50)) {
        if (mapToSameContig(clt->flag))
            return;
        fprintf(stderr, "\n[ Insertion region of %d_%d was refined to ===> CLT_LEFT_FLANK_MAP ]\tPrimary definition(flag tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\t%d ===> Secondary definition(tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\n", clt->tid, clt->idx, clt->flag, clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, tid1, leftMost, tid2, rightMost);
        clt->flag &= 0xffffffc3;
        clt->flag |= CLT_LEFT_FLANK_MAP;
        goto END;
    }

    // Only right end can be refined
    if ((tid1 < 0 || distToBp1 > 50) && (tid2 >= 0 && distToBp2 <= 50)) {
        if (mapToSameContig(clt->flag))
            return;
        fprintf(stderr, "\n[ Insertion region of %d_%d was refined to ===> CLT_RIGHT_FLANK_MAP ]\tPrimary definition(flag tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\t%d ===> Secondary definition(tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\n", clt->tid, clt->idx, clt->flag, clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, tid1, leftMost, tid2, rightMost);
        clt->flag &= 0xffffffc3;
        clt->flag |= CLT_RIGHT_FLANK_MAP;
        goto END;
    }

    // Both ends can be refiend
    if ((tid1 >= 0 && tid2 >= 0) && (distToBp1 <= 50 && distToBp2 <= 50)) {
        if (leftMost > rightMost - 100)
            return;
        fprintf(stderr, "\n[ Insertion region of %d_%d was refined to ===> CLT_SAME_FLANK_MAP ]\tPrimary definition(flag tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\t%d ===> Secondary definition(tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\n", clt->tid, clt->idx, clt->flag, clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, tid1, leftMost, tid2, rightMost);
        clt->flag &= 0xffffffc3;
        clt->flag |= (tid1 == tid2) ? CLT_SAME_FLANK_MAP : CLT_DIFF_FLANK_MAP;
        goto END;
    }

    if (isFlankMapped(clt->flag)) {
        fprintf(stderr, "\n[ Failed to refine insertion region of %d_%d, but primary definition succeeded ]\tPrimary definition(flag tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\t%d <=== Secondary definition(tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\n", clt->tid, clt->idx, clt->flag, clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, tid1, leftMost, tid2, rightMost);
        return;
    } else {
        fprintf(stderr, "\n[ Failed to define insertion region of %d_%d, check tmp_anno/%d_%d_FlankToAssm.bam to verify primary definition, check tmp_anno/%d_%d_AssmToRef.bam to verify secondary definition ]\tPrimary definition(flag tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\t%d ==== Secondary definition(tid1 leftMost tid2 rightMost): %d\t%d\t%d\t%d\n", clt->tid, clt->idx, clt->tid, clt->idx, clt->tid, clt->idx, clt->flag, clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, tid1, leftMost, tid2, rightMost);
    }

    END:
        clt->tid1 = tid1;
        clt->leftMost = leftMost;
        clt->tid2 = tid2;
        clt->rightMost = rightMost;
}

/// @brief output insertion-seq for annotation
void outputInsSeq(faidx_t *assmFa, Cluster *clt)
{
    char *insSeq = getInsSeq(assmFa, clt);
    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);

    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">%d_%d_%d_%d\n%s\n", clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, insSeq);
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

// /// @brief output flank-seqs of insertion-seq from assembled-contig for re-defining insertion region
// void outputAssmFlank(faidx_t *assmFa, Cluster *clt)
// {
//     int insLen = isLeftFlankMapped(clt->flag) ? (faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1)) - clt->leftMost) : clt->rightMost;
//     if (insLen < 250) {
//         outputInsSeq(assmFa, clt);
//         return;
//     }

//     char *leftSeq = NULL, *rightSeq = NULL;
//     hts_pos_t seqLen;

//     if (isLeftFlankMapped(clt->flag)) {
//         int assmLen = faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1));
//         leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost-200, clt->leftMost-1, &seqLen);
//         rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), assmLen-200, assmLen-1, &seqLen);
//     }

//     if (isRightFlankMapped(clt->flag)) {
//         leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, 199, &seqLen);
//         rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, clt->rightMost+199, &seqLen);
//     }

//     char outFn[100] = {'\0'};
//     sprintf(outFn, "tmp_anno/%d_%d_assmFlank.fa", clt->tid, clt->idx);
//     FILE *fp = fopen(outFn, "w");
//     fprintf(fp, ">0\n%s\n", leftSeq);
//     fprintf(fp, ">1\n%s\n", rightSeq);
//     fclose(fp);

//     if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
//     if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
// }

// /// @brief refine and ouput insertion-seq for single-flank-mapped cases
// void reExtractIns(Cluster *clt)
// {
//     char assmFn[100] = {'\0'};
//     sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
//     faidx_t *assmFa = fai_load((const char *)assmFn);
//     reSetInsRegion(clt, assmFa);
//     outputInsSeq(assmFa, clt);

//     if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
//     return;
// }

// /// @brief re-set insertion-seq region
// void reSetInsRegion(Cluster *clt, faidx_t *assmFa)
// {
//     char inputFn[100] = {'\0'};
//     sprintf(inputFn, "tmp_anno/%d_%d_AssmFlankToLocal.bam", clt->tid, clt->idx);
//     htsFile *inputBam = sam_open(inputFn, "rb");
//     sam_hdr_t *header = sam_hdr_read(inputBam);
//     bam1_t *bam = bam_init1();
//     int leftEnd = 0, rightStart = 0, leftLen = 0, rightLen = 0;
//     int localStart = atoi(sam_hdr_tid2name(header, 0));
//     uint32_t leftCigar1 = 0, leftCigar2 = 0;
//     uint32_t rightCigar1 = 0, rightCigar2 = 0;

//     while (1)
//     {
//         int retValue = bam_read1(inputBam->fp.bgzf, bam);
//         if (retValue < 0)
//             break;
//         if (bamIsInvalid(bam) || bamIsSup(bam) || bam_is_rev(bam))
//             continue;

//         int numCigar = bam->core.n_cigar;
//         uint32_t *cigarArr = bam_get_cigar(bam);

//         if (isLeftFlank(bam)) {
//             leftEnd = bam_endpos(bam);
//             leftLen = bam->core.l_qseq;
//             leftCigar1 = cigarArr[0];
//             leftCigar2 = cigarArr[numCigar - 1];
//         } else {
//             rightStart = bam->core.pos;
//             rightLen = bam->core.l_qseq;
//             rightCigar1 = cigarArr[0];
//             rightCigar2 = cigarArr[numCigar - 1];
//         }
//     }
    
//     if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
//     if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
//     if (header != NULL) {sam_hdr_destroy(header); header=NULL;}

//     if (leftEnd == 0 || rightStart == 0 )
//         return;

//     if (isClipInFlank(leftCigar1, 50) || isClipInFlank(rightCigar2, 50))
//         return;

//     int leftDelta = MIN(abs(leftEnd + localStart - clt->refStart), abs(leftEnd + localStart - clt->refEnd));
//     int rightDelta = MIN(abs(rightStart + localStart - clt->refStart), abs(rightStart + localStart - clt->refEnd));
//     if (leftDelta > 50 || rightDelta > 50)
//         return;

//     if (abs(rightStart - leftEnd) > 50)
//         return;

//     if (isLeftFlankMapped(clt->flag)) {
//         clt->tid2 = clt->tid1;
//         clt->rightMost = faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1));
//         clt->rightMost -= (rightLen - bam_cigar_oplen(rightCigar1));
//     }

//     if (isRightFlankMapped(clt->flag)) {
//         clt->tid1 = clt->tid2;
//         clt->leftMost = leftLen - bam_cigar_oplen(leftCigar2);
//     }

//     clt->flag &= 0xfffffff3;
//     clt->flag |= CLT_SAME_FLANK_MAP;
// }

/*******************
 *** Cluster I/O ***
 *******************/

/// @brief Output formated cluster records
void outputClt(Cluster *cltArr, int startIdx, int endIdx, const char *refFn, const char *teFn)
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
        Cluster *clt = &cltArr[i];
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
    if (hasTSD(clt->flag))
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

    if (isLeftFlankMapped(clt->flag)) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        **rightSeq = '.';
        goto END;
    }

    if (isRightFlankMapped(clt->flag)) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
        **leftSeq = '.';
        goto END;
    }
    
    *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
    *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
    goto END;

    END:
    *insSeq = getAdjustSeq(assmFa, clt);
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}

/// @brief get adjusted insertion-seq
char *getAdjustSeq(faidx_t *assmFa, Cluster *clt)
{
    char *insSeq = getInsSeq(assmFa, clt);
    if (isBothFlankMapped(clt->flag))
        return insSeq;

    int insLen = strlen(insSeq);
    int adjustLen = strlen(insSeq) + 50;
    char *adjustSeq = (char *)malloc((adjustLen + 1) * sizeof(char));
    char *temp = adjustSeq;

    if (isLeftFlankMapped(clt->flag)) {
        memcpy(temp, insSeq, insLen); temp += insLen;
        memset(temp, 'N', 50); temp += 50;
        adjustSeq[adjustLen] = '\0';

        if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
        return adjustSeq;
    }

    if (isRightFlankMapped(clt->flag)) {
        memset(temp, 'N', 50); temp += 50;
        memcpy(temp, insSeq, insLen); temp += insLen;
        adjustSeq[adjustLen] = '\0';

        if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
        return adjustSeq;
    }

    return NULL;
}