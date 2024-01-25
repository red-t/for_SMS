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
void getTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize)
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

/// @brief define flank region on ref-genome
int getFlankRegion(Cluster *cluster, int *leftStart, int *leftEnd, int *rightStart, int *rightEnd)
{
    if (cluster->refEnd - cluster->refStart == 1)
    {
        *leftStart = cluster->refEnd - 450;
        *leftEnd = cluster->refEnd + 50;
        *rightStart = cluster->refStart - 50;
        *rightEnd = cluster->refStart + 450;
        return 0;
    }

    *leftStart = cluster->refEnd - 500;
    *leftEnd = cluster->refEnd;
    *rightStart = cluster->refStart;
    *rightEnd = cluster->refStart + 500;
    return 0;
}

/// @brief extract and output flank sequence for all clusters
int outputRefFlankSeqs(char *refFn, Cluster *cltArray, int startIdx, int endIdx)
{
    hts_pos_t seqLen;
    int leftStart, leftEnd, rightStart, rightEnd;
    faidx_t *refFa = fai_load((const char*)refFn);
    char *outFn = (char *)malloc(100 * sizeof(char));

    for (int i = startIdx; i < endIdx; i++)
    {
        Cluster *cluster = &cltArray[i];
        getFlankRegion(cluster, &leftStart, &leftEnd, &rightStart, &rightEnd);
        const char *chrom = faidx_iseq(refFa, cluster->tid);
        char *leftSeq = faidx_fetch_seq64(refFa, chrom, leftStart, leftEnd, &seqLen);
        char *rightSeq = faidx_fetch_seq64(refFa, chrom, rightStart, rightEnd, &seqLen);

        sprintf(outFn, "tmp_anno/%d_%d_flank.fa", cluster->tid, cluster->idx);
        FILE *fp = fopen(outFn, "w");
        fprintf(fp, ">0\n%s\n", leftSeq);
        fprintf(fp, ">1\n%s\n", rightSeq);
        fclose(fp);
    }

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
    return 0;
}

/*****************************
 *** Insertion Sequence IO ***
 *****************************/

/// @brief adjust insertion sequence region
int adjustInsRegion(Cluster *cluster, int leftTid, int rightTid, uint32_t leftCigar, uint32_t rightCigar, int leftLen)
{
    // Only right-flank mapped
    if (leftTid < 0 && rightTid >= 0) {
        cluster->flag |= CLT_SINGLE_FLANK_MAP;
        return rightTid;
    }

    // Only left-flank mapped
    if (leftTid >= 0 && rightTid < 0) {
        cluster->insEnd = leftLen;
        cluster->flag |= CLT_SINGLE_FLANK_MAP;
        return leftTid;
    }

    // Two flanks map to different contig
    if (leftTid != rightTid)
        return 0;

    // Two flanks map to same contig
    if (cluster->insStart < cluster->insEnd) {
        cluster->flag |= CLT_BOTH_FLANK_MAP;
        return leftTid;
    }

    // Two flanks map to same contig, but only right is qualified
    if (bam_cigar_op(leftCigar) == BAM_CSOFT_CLIP && bam_cigar_oplen(leftCigar) > 400) {
        cluster->insStart = 0;
        cluster->flag |= CLT_SINGLE_FLANK_MAP;
        return rightTid;
    }

    // Two flanks map to same contig, but only left is qualified
    if (bam_cigar_op(rightCigar) == BAM_CSOFT_CLIP && bam_cigar_oplen(rightCigar) > 400) {
        cluster->insEnd = leftLen;
        cluster->flag |= CLT_SINGLE_FLANK_MAP;
        return leftTid;
    }

    return 0;
}

/// @brief define insertion sequence region
int getInsRegion(Cluster *cluster)
{
    char *inputFn = (char *)malloc(100 * sizeof(char));
    sprintf(inputFn, "tmp_anno/%d_%d_FlankToAssm.bam", cluster->tid, cluster->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bamRecord = bam_init1();
    int leftTid = -1, rightTid = -2, leftLen = 0;
    uint32_t leftCigar = 0, rightCigar = 0;

    while (1)
    {
        int returnValue = bam_read1(inputBam->fp.bgzf, bamRecord);
        if (returnValue < 0)
            break;

        if (bamIsInvalid(bamRecord) || bamIsSup(bamRecord) || bam_is_rev(bamRecord))
            continue;

        int numCigar = bamRecord->core.n_cigar;
        uint32_t *cigarArray = bam_get_cigar(bamRecord);
        if (isLeftFlank(bamRecord)) {
            cluster->insStart = bam_endpos(bamRecord);
            leftTid = bamRecord->core.tid;
            leftLen = sam_hdr_tid2len(header, leftTid);
            leftCigar = cigarArray[0];
        } else {
            cluster->insEnd = bamRecord->core.pos;
            rightTid = bamRecord->core.tid;
            rightCigar = cigarArray[numCigar-1];
        }
    }

    if (bamRecord != NULL) {bam_destroy1(bamRecord); bamRecord=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}

    return adjustInsRegion(cluster, leftTid, rightTid, leftCigar, rightCigar, leftLen);
}

/// @brief extract and output insertion sequence for single cluster
int outputInsSeq(Cluster *cluster)
{
    int contigTid = getInsRegion(cluster);
    if ((cluster->flag & CLT_SINGLE_FLANK_MAP) == 0 && (cluster->flag & CLT_BOTH_FLANK_MAP) == 0)
        return 0;

    hts_pos_t seqLen;
    char *refFn = (char *)malloc(100 * sizeof(char));
    sprintf(refFn, "tmp_assm/tmp.%d_%d_assembled.fa", cluster->tid, cluster->idx);
    faidx_t *refFa = fai_load((const char *)refFn);
    const char *contig = faidx_iseq(refFa, contigTid);
    char *insSeq = faidx_fetch_seq64(refFa, contig, cluster->insStart, cluster->insEnd, &seqLen);

    char *outFn = (char *)malloc(100 * sizeof(char));
    sprintf(outFn, "tmp_anno/%d_%d_insertion.fa", cluster->tid, cluster->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">%d-%d\n%s\n", cluster->insStart, cluster->insEnd, insSeq);
    fclose(fp);

    if (refFn != NULL) {free(refFn); refFn=NULL;}
    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
    if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
    if (outFn != NULL) {free(outFn); outFn=NULL;}
    return 0;
}