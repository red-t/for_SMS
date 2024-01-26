#include "anno_utils.h"

/// @brief get all annotate records by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx)
{
    char *inputFn = (char *)malloc(100 * sizeof(char));
    sprintf(inputFn, "tmp_anno/%d_%d_InsToTE.bam", cluster->tid, cluster->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bamRecord = bam_init1();

    int numAnno = 0;
    int leftMost = INT_MAX, rightMost = 0;
    int leftIdx = 0, rightIdx = 0;
    while (1)
    {
        int returnValue = bam_read1(inputBam->fp.bgzf, bamRecord);
        if (returnValue < 0)
            break;
        if (bamIsInvalid(bamRecord))
            continue;

        initAnno(bamRecord, &annoArray[numAnno], idx);
        if (annoArray[numAnno].queryStart < leftMost) {
            leftMost = annoArray[numAnno].queryStart;
            leftIdx = numAnno;
        }
        if (annoArray[numAnno].queryEnd > rightMost) {
            rightMost = annoArray[numAnno].queryEnd;
            rightIdx = numAnno;
        }
        numAnno++;
    }

    if (numAnno > 0)
        cluster->flag |= CLT_TE_MAP;

    int prevNumAnno = numAnno;
    numAnno = annoPolyA(cluster, idx, annoArray, numAnno, leftMost, rightMost, leftIdx, rightIdx);
    if (numAnno - prevNumAnno > 0)
        cluster->flag |= CLT_POLYA;

    if (bamRecord != NULL) {bam_destroy1(bamRecord); bamRecord=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}
    return numAnno;
}

/// @brief init single annotate record from bamRecord
void initAnno(bam1_t *bamRecord, Anno *anno, int idx)
{
    int numCigar = bamRecord->core.n_cigar;
    uint32_t *cigarArray = bam_get_cigar(bamRecord);

    int queryStart = 0;
    if (bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP)
        queryStart = bam_cigar_oplen(cigarArray[0]);

    int queryEnd = bamRecord->core.l_qseq;
    if (bam_cigar_op(cigarArray[numCigar - 1]) == BAM_CSOFT_CLIP)
        queryEnd -= bam_cigar_oplen(cigarArray[numCigar - 1]);

    int i, refLen;
    for (i = refLen = 0; i < numCigar; i++)
        if (bam_cigar_type(bam_cigar_op(cigarArray[i])) & 2)
            refLen += bam_cigar_oplen(cigarArray[i]);

    anno->idx = idx;
    anno->queryStart = queryStart;
    anno->queryEnd = queryEnd;
    anno->strand = bam_is_rev(bamRecord);
    anno->tid = bamRecord->core.tid;
    anno->refStart = bamRecord->core.pos;
    anno->refEnd = bamRecord->core.pos + refLen;

    if (bam_is_rev(bamRecord))
    {
        anno->queryStart = bamRecord->core.l_qseq - queryEnd;
        anno->queryEnd = bamRecord->core.l_qseq - queryStart;
    }
}

/// @brief annotate polyA/polyT
int annoPolyA(Cluster *cluster, int idx, Anno *annoArray, int numAnno, int leftMost, int rightMost, int leftIdx, int rightIdx)
{
    char *insFn = (char *)malloc(100 * sizeof(char));
    sprintf(insFn, "tmp_anno/%d_%d_insertion.fa", cluster->tid, cluster->idx);
    faidx_t *insFa = fai_load((const char *)insFn);
    const char *insID = faidx_iseq(insFa, 0);
    int insLen = faidx_seq_len64(insFa, insID);
    char *flankSeq = NULL;
    hts_pos_t seqLen;

    if (leftMost > 10) {
        flankSeq = faidx_fetch_seq64(insFa, insID, 0, leftMost, &seqLen);
        numAnno = getPolyA(flankSeq, seqLen, idx, 0, annoArray, numAnno, rightMost, leftIdx, rightIdx);
    }

    if ((insLen - rightMost) > 10) {
        flankSeq = faidx_fetch_seq64(insFa, insID, rightMost, insLen, &seqLen);
        numAnno = getPolyA(flankSeq, seqLen, idx, 1, annoArray, numAnno, rightMost, leftIdx, rightIdx);
    }

    if (insFn != NULL) {free(insFn); insFn=NULL;}
    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    if (flankSeq != NULL) {free(flankSeq); flankSeq=NULL;}
    return numAnno;
}

/// @brief find polyA/polyT region and init single annotate records
int getPolyA(char *flankSeq, int seqLen, int idx, int isA, Anno *annoArray, int numAnno, int rightMost, int leftIdx, int rightIdx)
{
    uint8_t targetChar = (isA) ? 65 : 84; // A OR T
    int thisLen = 0, maxLen = 0;
    int thisSum = 0, maxSum = 0;
    int thisNum = 0, maxNum = 0;
    int end = 0;

    for (int i = 0; i < seqLen; i++)
    {
        thisLen++;
        // Is a/A OR t/T
        if (((flankSeq[i] | 0x20) & 0x5f) == targetChar) {
            thisSum++;
            thisNum++;
        }
        else
            thisSum--;
        
        if (thisSum > maxSum) {
            maxLen = thisLen;
            maxSum = thisSum;
            maxNum = thisNum;
            end = i + 1;
        }

        if (thisSum < 0) {
            thisSum = 0;
            thisNum = 0;
            thisLen = 0;
        }
    }

    // polyA at the end, but right-most TE is reverse
    if (isA && annoArray[rightIdx].strand)
        return numAnno;

    // polyT at the start, but left-most TE is forward
    if ((!isA) && (!annoArray[rightIdx].strand))
        return numAnno;

    if(maxLen < 10)
        return numAnno;

    float fracA = (float)maxNum / maxLen;
    if (fracA < 0.8)
        return numAnno;

    end = (isA) ? (end + rightMost) : end;
    annoArray[numAnno].idx = idx;
    annoArray[numAnno].queryStart = end - maxLen;
    annoArray[numAnno].queryEnd = end;
    annoArray[numAnno].strand = (isA) ? 0 : 1;
    annoArray[numAnno].tid = (isA) ? -1 : -2;
    numAnno++;
    return numAnno;
}
