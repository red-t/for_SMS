#include "anno_utils.h"

/// @brief Initiate PolyA
PolyA initPolyA(int idx, int cltTid, int cltIdx)
{
    PolyA polyA;
    polyA.leftMost = INT_MAX;
    polyA.rightMost = 0;
    polyA.idx = idx;
    polyA.cltTid = cltTid;
    polyA.cltIdx = cltIdx;
    return polyA;
}

/// @brief Record single TE annotation
void initAnno(bam1_t *bam, Anno *anno, int idx, int cltTid, int cltIdx)
{
    int numCigar = bam->core.n_cigar;
    uint32_t *cigarArray = bam_get_cigar(bam);

    int queryStart = 0;
    if (bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP)
        queryStart = bam_cigar_oplen(cigarArray[0]);

    int queryEnd = bam->core.l_qseq;
    if (bam_cigar_op(cigarArray[numCigar - 1]) == BAM_CSOFT_CLIP)
        queryEnd -= bam_cigar_oplen(cigarArray[numCigar - 1]);

    int i, refLen;
    for (i = refLen = 0; i < numCigar; i++)
        if (bam_cigar_type(bam_cigar_op(cigarArray[i])) & 2)
            refLen += bam_cigar_oplen(cigarArray[i]);

    anno->idx = idx;
    anno->cltTid = cltTid;
    anno->cltIdx = cltIdx;
    anno->queryStart = queryStart;
    anno->queryEnd = queryEnd;
    anno->strand = bam_is_rev(bam);
    anno->tid = bam->core.tid;
    anno->refStart = bam->core.pos;
    anno->refEnd = bam->core.pos + refLen;

    if (bam_is_rev(bam))
    {
        anno->queryStart = bam->core.l_qseq - queryEnd;
        anno->queryEnd = bam->core.l_qseq - queryStart;
    }
}

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx)
{
    char *inputFn = (char *)malloc(100 * sizeof(char));
    sprintf(inputFn, "tmp_anno/%d_%d_InsToTE.bam", cluster->tid, cluster->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int numAnno = 0;
    PolyA polyA = initPolyA(idx, cluster->tid, cluster->idx);
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam))
            continue;

        initAnno(bam, &annoArray[numAnno], idx, cluster->tid, cluster->idx);
        if (annoArray[numAnno].queryStart < polyA.leftMost) {
            polyA.leftMost = annoArray[numAnno].queryStart;
            polyA.leftIdx = numAnno;
        }
        if (annoArray[numAnno].queryEnd > polyA.rightMost) {
            polyA.rightMost = annoArray[numAnno].queryEnd;
            polyA.rightIdx = numAnno;
        }
        numAnno++;
    }

    if (numAnno > 0)
        cluster->flag |= CLT_TE_MAP;

    int numTE = numAnno;
    numAnno = annoPolyA(cluster, annoArray, numAnno, polyA);
    if (numAnno - numTE > 0)
        cluster->flag |= CLT_POLYA;

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}
    return numAnno;
}

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *cluster, Anno *annoArray, int numAnno, PolyA polyA)
{
    char *insFn = (char *)malloc(100 * sizeof(char));
    sprintf(insFn, "tmp_anno/%d_%d_insertion.fa", cluster->tid, cluster->idx);
    faidx_t *insFa = fai_load((const char *)insFn);
    const char *insID = faidx_iseq(insFa, 0);
    char *flankSeq = NULL;
    hts_pos_t seqLen;

    if (polyA.leftMost > 10) {
        flankSeq = faidx_fetch_seq64(insFa, insID, 0, polyA.leftMost, &seqLen);
        polyA.isA = 0;
        polyA.seqLen = seqLen;
        numAnno = getPolyA(flankSeq, annoArray, numAnno, polyA);
    }

    int insLen = faidx_seq_len64(insFa, insID);
    if ((insLen - polyA.rightMost) > 10) {
        flankSeq = faidx_fetch_seq64(insFa, insID, polyA.rightMost, insLen, &seqLen);
        polyA.isA = 1;
        polyA.seqLen = seqLen;
        numAnno = getPolyA(flankSeq, annoArray, numAnno, polyA);
    }

    if (insFn != NULL) {free(insFn); insFn=NULL;}
    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    if (flankSeq != NULL) {free(flankSeq); flankSeq=NULL;}
    return numAnno;
}

/// @brief Find and record single polyA/polyT
int getPolyA(char *flankSeq, Anno *annoArray, int numAnno, PolyA polyA)
{
    uint8_t targetChar = (polyA.isA) ? 65 : 84; // A OR T
    int thisLen = 0, maxLen = 0;
    int thisSum = 0, maxSum = 0;
    int thisNum = 0, maxNum = 0;
    int end = 0;

    for (int i = 0; i < polyA.seqLen; i++)
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
    if (polyA.isA && annoArray[polyA.rightIdx].strand)
        return numAnno;

    // polyT at the start, but left-most TE is forward
    if ((!polyA.isA) && (!annoArray[polyA.leftIdx].strand))
        return numAnno;

    if(maxLen < 10)
        return numAnno;

    float fracA = (float)maxNum / maxLen;
    if (fracA < 0.8)
        return numAnno;

    end = (polyA.isA) ? (end + polyA.rightMost) : end;
    annoArray[numAnno].idx = polyA.idx;
    annoArray[numAnno].cltTid = polyA.cltTid;
    annoArray[numAnno].cltIdx = polyA.cltIdx;
    annoArray[numAnno].queryStart = end - maxLen;
    annoArray[numAnno].queryEnd = end;
    annoArray[numAnno].strand = (polyA.isA) ? 0 : 1;
    annoArray[numAnno].tid = (polyA.isA) ? -1 : -2;
    numAnno++;
    return numAnno;
}

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annoTsd(Cluster *cluster)
{
    if (!isBothFlankMapped(cluster->flag))
        return;

    char *inputFn = (char *)malloc(100 * sizeof(char));
    sprintf(inputFn, "tmp_anno/%d_%d_TsdToLocal.bam", cluster->tid, cluster->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int leftEnd = -1, rightStart = -1;
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam) || bamIsSup(bam) || bam_is_rev(bam))
            continue;

        if (isLeftFlank(bam))
            leftEnd = bam_endpos(bam);
        else
            rightStart = bam->core.pos;
    }

    setTsd(cluster, atoi(sam_hdr_tid2name(header, 0)), leftEnd, rightStart);

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}
}

/// @brief Find TSD and refine breakpoint
void setTsd(Cluster *cluster, int localStart, int leftEnd, int rightStart)
{
    if (leftEnd < 0 && rightStart < 0)
        return;
    
    if (leftEnd < 0) {
        cluster->refEnd = localStart + rightStart;
        cluster->refStart = cluster->refEnd - 1;
        return;
    }

    if (rightStart < 0) {
        cluster->refStart = localStart + leftEnd;
        cluster->refEnd = cluster->refStart + 1;
        return;
    }

    if (rightStart < leftEnd && (leftEnd - rightStart) < 50) {
        cluster->flag |= CLT_TSD;
        cluster->refStart = localStart + rightStart;
        cluster->refEnd = cluster->refStart + 1;
        cluster->tsdStart = localStart + rightStart;
        cluster->tsdEnd = localStart + leftEnd;
        return;
    }

    if (rightStart > leftEnd) {
        cluster->refStart = localStart + leftEnd;
        cluster->refEnd = localStart + rightStart;
    }
}

/// @brief Output annotation records
void outPutAnno(Anno *annoArray, int numAnno, const char *outFn)
{
    char *queryTmp = malloc(100 * sizeof(char));
    char *refTmp = malloc(100 * sizeof(char));
    char *queryStr = malloc(500 * sizeof(char));
    char *refStr = malloc(500 * sizeof(char));
    memset(queryStr, '\0', 500);
    memset(refStr, '\0', 500);
    
    int prevIdx = annoArray[0].idx;
    FILE *fp = fopen(outFn, "w");
    for (int i = 0; i < numAnno; i++)
    {
        if (annoArray[i].idx != prevIdx) {
            queryStr[strlen(queryStr)-1] = '\0';
            refStr[strlen(refStr)-1] = '\0';
            fprintf(fp, "%d-%d\t%s\t%s\n", annoArray[i-1].cltTid, annoArray[i-1].cltIdx, queryStr, refStr);
            prevIdx = annoArray[i].idx;
            memset(queryStr, '\0', strlen(queryStr));
            memset(refStr, '\0', strlen(refStr));
        }

        char strand = (annoArray[i].strand == 0) ? '+' : '-';
        sprintf(queryTmp, "%c:%d-%d,", strand, annoArray[i].queryStart, annoArray[i].queryEnd);
        sprintf(refTmp, "%d:%d-%d,", annoArray[i].tid, annoArray[i].refStart, annoArray[i].refEnd);
        strcat(queryStr, queryTmp);
        strcat(refStr, refTmp);
    }
    // Output final anno record
    queryStr[strlen(queryStr)-1] = '\0';
    refStr[strlen(refStr)-1] = '\0';
    fprintf(fp, "%d-%d\t%s\t%s\n", annoArray[numAnno-1].cltTid, annoArray[numAnno-1].cltIdx, queryStr, refStr);
    fclose(fp);

    if (queryStr != NULL) {free(queryStr); queryStr=NULL;}
    if (refStr != NULL) {free(refStr); refStr=NULL;}
    if (queryTmp != NULL) {free(queryTmp); queryTmp=NULL;}
    if (refTmp != NULL) {free(refTmp); refTmp=NULL;}
}