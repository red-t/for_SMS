#include <stdlib.h>
#include "anno_utils.h"

/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/

/// @brief Initiate PolyA
PolyA initPolyA(int idx)
{
    PolyA polyA;
    polyA.leftMost = INT_MAX;
    polyA.rightMost = 0;
    polyA.idx = idx;
    return polyA;
}

/// @brief Record single TE annotation
void initAnno(bam1_t *bam, sam_hdr_t *header, Cluster *clt, Anno *anno, int idx)
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
    anno->cltTid = clt->tid;
    anno->cltIdx = clt->idx;
    anno->queryStart = queryStart;
    anno->queryEnd = queryEnd;
    anno->strand = bam_is_rev(bam);
    anno->tid = bam->core.tid;
    anno->refStart = bam->core.pos;
    anno->refEnd = bam->core.pos + refLen;

    if (bam_is_rev(bam)) {
        anno->queryStart = bam->core.l_qseq - queryEnd;
        anno->queryEnd = bam->core.l_qseq - queryStart;
    }

    int teLen = sam_hdr_tid2len(header, anno->tid);
    int truncSize = (((float)teLen * 0.08) < 100) ? (teLen * 0.08) : 100;
    if (anno->refStart < truncSize)
        anno->flag |= CLT_5P_FULL;
    if ((teLen - anno->refEnd) < truncSize)
        anno->flag |= CLT_3P_FULL;
    
    if (anno->tid == clt->repTid)
        clt->flag |= CLT_SELF_TO_SELF;
}

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *clt, Anno *annoArray, int idx)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_InsToTE.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int numAnno = 0;
    PolyA polyA = initPolyA(idx);
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam))
            continue;

        initAnno(bam, header, clt, &annoArray[numAnno], idx);
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

    if (numAnno == 0)
        goto END;

    if (numAnno > 1)
        clt->flag |= CLT_MULTI_TE;
    clt->flag |= CLT_TE_MAP;
    numAnno = annoPolyA(clt, annoArray, numAnno, &polyA);
    outputTsdSeq(clt, &polyA, annoArray, numAnno);
    goto END;

    END:
    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    return numAnno;
}

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *clt, Anno *annoArray, int numAnno, PolyA *polyA)
{
    char insFn[100] = {'\0'};
    sprintf(insFn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);
    faidx_t *insFa = fai_load((const char *)insFn);
    const char *insID = faidx_iseq(insFa, 0);
    char *flankSeq = NULL;
    hts_pos_t seqLen;

    if (polyA->leftMost >= 5) {
        flankSeq = faidx_fetch_seq64(insFa, insID, 0, polyA->leftMost-1, &seqLen);
        polyA->isA = 0;
        polyA->seqLen = seqLen;
        numAnno = setPolyA(flankSeq, annoArray, clt, numAnno, polyA);
    }

    clt->insLen = faidx_seq_len64(insFa, insID);
    if ((clt->insLen - polyA->rightMost) >= 5) {
        flankSeq = faidx_fetch_seq64(insFa, insID, polyA->rightMost, clt->insLen-1, &seqLen);
        polyA->isA = 1;
        polyA->seqLen = seqLen;
        numAnno = setPolyA(flankSeq, annoArray, clt, numAnno, polyA);
    }

    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    if (flankSeq != NULL) {free(flankSeq); flankSeq=NULL;}
    return numAnno;
}

/// @brief Find and record single polyA/polyT
int setPolyA(char *flankSeq, Anno *annoArray, Cluster *clt, int numAnno, PolyA *polyA)
{
    // polyA at right, but right-most TE is reverse
    if (polyA->isA && annoArray[polyA->rightIdx].strand)
        return numAnno;
    // polyT at left, but left-most TE is forward
    if (!polyA->isA && !annoArray[polyA->leftIdx].strand)
        return numAnno;

    int thisLen = 0, maxLen = 0;
    int thisSum = 0, maxSum = 0;
    int thisNum = 0, maxNum = 0;
    int numOther = 0, stop = 0;

    uint8_t targetChar = (polyA->isA) ? 65 : 84;
    int start = (polyA->isA) ? 0 : polyA->seqLen-1;
    int end = (polyA->isA) ? polyA->seqLen : -1;
    int step = (polyA->isA) ? 1 : -1;

    for (int i = start; i != end; i += step)
    {
        thisLen++;
        // Is a/A OR t/T
        if (((flankSeq[i] | 0x20) & 0x5f) == targetChar) {
            thisSum++;
            thisNum++;
        } else {
            thisSum--;
            numOther++;
        }
        if (thisSum > maxSum) {
            maxLen = thisLen;
            maxSum = thisSum;
            maxNum = thisNum;
            stop = i;
        }
        if (thisSum < 0 || numOther > 3)
            thisSum = thisNum = thisLen = numOther = 0;
    }

    if(maxLen < 5)
        return numAnno;
    if (((float)maxNum / maxLen) < 0.8)
        return numAnno;

    if (polyA->isA) {
        annoArray[numAnno].queryStart = stop + polyA->rightMost + 1 - maxLen;
        annoArray[numAnno].queryEnd = stop + polyA->rightMost + 1;
        annoArray[numAnno].strand = 0;
        annoArray[numAnno].tid = -1;
        polyA->rightMost = annoArray[numAnno].queryEnd;
    } else {
        annoArray[numAnno].queryStart = stop;
        annoArray[numAnno].queryEnd = stop + maxLen;
        annoArray[numAnno].strand = 1;
        annoArray[numAnno].tid = -2;
        polyA->leftMost = annoArray[numAnno].queryStart;
    }
    annoArray[numAnno].idx = polyA->idx;
    annoArray[numAnno].cltTid = clt->tid;
    annoArray[numAnno].cltIdx = clt->idx;

    return ++numAnno;
}

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *clt, PolyA *polyA, Anno *annoArray, int numAnno)
{
    if (!isBothFlankMapped(clt->flag))
        return;

    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);
    hts_pos_t seqLen;

    int leftStart = clt->leftMost - 100;
    int leftEnd = clt->leftMost + polyA->leftMost;
    int rightStart = clt->rightMost - (clt->insLen - polyA->rightMost);
    int rightEnd = clt->rightMost + 100;
    char *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), leftStart, leftEnd-1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), rightStart, rightEnd-1, &seqLen);

    if ((leftEnd > faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1))) || (rightStart < 0))
        goto END;

    adjustAnno(annoArray, numAnno, -polyA->leftMost);
    clt->leftMost = leftEnd;
    clt->rightMost = rightStart;
    clt->insLen = clt->rightMost - clt->leftMost;

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_tsd.fa", clt->tid, clt->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);
    goto END;

    END:
    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}

/// @brief Adjust annotation position
void adjustAnno(Anno *annoArray, int numAnno, int leftDelta)
{
    for (int i = 0; i < numAnno; i++)
    {
        annoArray[i].queryStart += leftDelta;
        annoArray[i].queryEnd += leftDelta;
    }
}

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annoTsd(Cluster *clt, Anno *annoArray, int numAnno)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_TsdToLocal.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int leftEnd = -1, rightStart = -1;
    int leftDelta = 0, rightDelta = 0;
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam) || bam_is_rev(bam))
            continue;
        
        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArray = bam_get_cigar(bam);
        
        if (isLeftFlank(bam)) {
            if (bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP && bam_cigar_oplen(cigarArray[0]) > 90)
                continue;
            if (bam_cigar_op(cigarArray[numCigar-1]) == BAM_CSOFT_CLIP)
                leftDelta = bam_cigar_oplen(cigarArray[numCigar-1]);
            leftEnd = bam_endpos(bam);
        } else {
            if (bam_cigar_op(cigarArray[numCigar-1]) == BAM_CSOFT_CLIP && bam_cigar_oplen(cigarArray[numCigar-1]) > 90)
                continue;
            if (bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP)
                rightDelta = bam_cigar_oplen(cigarArray[0]);
            rightStart = bam->core.pos;
        }
    }

    setTsd(clt, atoi(sam_hdr_tid2name(header, 0)), leftEnd, rightStart);
    adjustAnno(annoArray, numAnno, leftDelta);
    clt->leftMost -= leftDelta;
    clt->rightMost += rightDelta;
    clt->insLen += leftDelta + rightDelta;

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
}

/// @brief Find TSD and refine breakpoint
void setTsd(Cluster *clt, int localStart, int leftEnd, int rightStart)
{
    if (leftEnd < 0 && rightStart < 0)
        return;
    
    if (leftEnd < 0) {
        clt->refEnd = localStart + rightStart;
        clt->refStart = clt->refEnd - 1;
        return;
    }

    if (rightStart < 0) {
        clt->refEnd = localStart + leftEnd;
        clt->refStart = clt->refEnd - 1;
        return;
    }

    if (rightStart < leftEnd) {
        clt->refStart = localStart + rightStart;
        clt->refEnd = localStart + leftEnd;
        clt->flag |= ((leftEnd - rightStart) < 50) ? CLT_TSD : 0;
        return;
    }
    
    clt->refStart = localStart + leftEnd - 1;
    clt->refEnd = localStart + rightStart;
}

/// @brief Set ins-seq structure based on annotations
void setInsStruc(Cluster *clt, Anno *annoArray, int numAnno, uint32_t *classArray)
{
    if (numAnno == 0)
        return;

    qsort(annoArray, numAnno, sizeof(Anno), compare);
    checkGap(annoArray, numAnno, clt);
    checkPolyA(annoArray, numAnno, clt);
    checkEnd(annoArray, numAnno, clt);
    checkTEClass(annoArray, numAnno, clt, classArray);
}

/// @brief Compare function for sorting annotations
int compare(const void *a, const void *b)
{
    Anno *pa = (Anno *)a;
    Anno *pb = (Anno *)b;
    if (pa->queryStart != pb->queryStart)
        return pa->queryStart - pb->queryStart;
    else
        return pa->queryEnd - pb->queryEnd;
}

/// @brief Check whether the ins-seq contains large gap
void checkGap(Anno *annoArray, int numAnno, Cluster *clt)
{
    int thisGap = 0, maxGap = annoArray[0].queryStart;
    for (int i = 1; i < numAnno; i++)
    {
        thisGap = annoArray[i].queryStart - annoArray[i-1].queryEnd;
        maxGap = (thisGap > maxGap) ? thisGap : maxGap;
    }
    thisGap = clt->insLen - annoArray[numAnno - 1].queryEnd;
    maxGap = (thisGap > maxGap) ? thisGap : maxGap;
    clt->flag |= (maxGap >= 1500) ? CLT_LARGE_GAP : 0;
}

/// @brief Check whether the ins-seq contains valid polyA
void checkPolyA(Anno *annoArray, int numAnno, Cluster *clt)
{
    if (numAnno == 1)
        return;

    int hasPolyT = (annoArray[0].tid == -2) ? 1 : 0;
    if (hasPolyT && is3PFull(annoArray[1].flag)) {
        int gap1 = annoArray[1].queryStart - annoArray[0].queryEnd;
        int gap2 = annoArray[0].queryStart;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }
    int hasPolyA = (annoArray[numAnno-1].tid == -1) ? 1 : 0;
    if (hasPolyA && is3PFull(annoArray[numAnno-2].flag)) {
        int gap1 = annoArray[numAnno-1].queryStart - annoArray[numAnno-2].queryEnd;
        int gap2 = clt->insLen - annoArray[numAnno-1].queryEnd;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }
}

/// @brief Check whether the ins-seq contains complete ends
void checkEnd(Anno *annoArray, int numAnno, Cluster *clt)
{
    int leftIdx = (annoArray[0].tid == -2) ? 1 : 0;
    int rightIdx = (annoArray[numAnno-1].tid == -1) ? numAnno-2 : numAnno-1;

    if (!isRevAnno(annoArray[leftIdx]) && is5PFull(annoArray[leftIdx].flag))
        clt->flag |= CLT_5P_FULL;
    if (isRevAnno(annoArray[leftIdx]) && is3PFull(annoArray[leftIdx].flag))
        clt->flag |= CLT_3P_FULL;

    if (!isRevAnno(annoArray[rightIdx]) && is3PFull(annoArray[rightIdx].flag))
        clt->flag |= CLT_3P_FULL;
    if (isRevAnno(annoArray[rightIdx]) && is5PFull(annoArray[rightIdx].flag))
        clt->flag |= CLT_5P_FULL;
}

/// @brief Check which TE class the insertion belongs to
void checkTEClass(Anno *annoArray, int numAnno, Cluster *clt, uint32_t *classArray)
{
    int leftIdx = (annoArray[0].tid == -2) ? 1 : 0;
    int rightIdx = (annoArray[numAnno-1].tid == -1) ? numAnno-2 : numAnno-1;

    int leftLen = annoArray[leftIdx].refEnd - annoArray[leftIdx].refStart;
    int rightLen = annoArray[rightIdx].refEnd - annoArray[rightIdx].refStart;
    int teTid = (leftLen > rightLen) ? annoArray[leftIdx].tid : annoArray[rightIdx].tid;

    clt->flag |= classArray[teTid];
}


/**********************
 *** Annotation I/O ***
 **********************/

/// @brief Output formated annotation records
void outputAnno(Anno *annoArray, int numAnno, int startIdx, const char *teFn)
{   
    faidx_t *teFa = fai_load(teFn);
    int numTe = faidx_nseq(teFa), strandFlag = 0, prevIdx = annoArray[0].idx;
    int *teTable = malloc(numTe * sizeof(int));
    memset(teTable, 0, numTe * sizeof(int));

    char queryTmp[100] = {'\0'}, refTmp[100] = {'\0'};
    char queryStr[1000] = {'\0'}, refStr[1000] = {'\0'};
    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_annoFormated.txt", startIdx);
    FILE *fp = fopen(outFn, "w");

    for (int i = 0; i < numAnno; i++)
    {
        if (annoArray[i].idx != prevIdx) {
            writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, annoArray[i-1]);
            prevIdx = annoArray[i].idx;
            strandFlag = 0;
            memset(teTable, 0, numTe * sizeof(int));
            memset(queryStr, '\0', 1000);
            memset(refStr, '\0', 1000);
        }
        formatSingleAnno(annoArray[i], queryTmp, refTmp, teFa, teTable, &strandFlag);
        strcat(queryStr, queryTmp);
        strcat(refStr, refTmp);
    }
    // Output final anno record
    writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, annoArray[numAnno-1]);
    fclose(fp);

    if (teFa != NULL) {fai_destroy(teFa); teFa = NULL;}
    if (teTable != NULL) {free(teTable); teTable = NULL;}
}

/// @brief Change single annotation record into specified format
void formatSingleAnno(Anno anno, char *queryTmp, char *refTmp, faidx_t *teFa, int *teTable, int *strandFlag)
{
    char strand = (anno.strand == 0) ? '+' : '-';
    sprintf(queryTmp, "%c:%d-%d,", strand, anno.queryStart, anno.queryEnd);

    if (anno.tid == -1)
        sprintf(refTmp, "PolyA:%d-%d,", anno.refStart, anno.refEnd);
    else if (anno.tid == -2)
        sprintf(refTmp, "PolyT:%d-%d,", anno.refStart, anno.refEnd);
    else
    {
        sprintf(refTmp, "%s:%d-%d,", faidx_iseq(teFa, anno.tid), anno.refStart, anno.refEnd);
        teTable[anno.tid] = 1;
        *strandFlag += (anno.strand == 0) ? 1 : -1;
    }
}

/// @brief Write annotation for a single cluster
void writeSingleCltAnno(int strandFlag, int numTe, int *teTable, faidx_t *teFa, char *queryStr, char *refStr, FILE *fp, Anno anno)
{
    char strand = getCltStrand(strandFlag);
    char *cltClass = getCltClass(numTe, teTable, teFa);
    queryStr[strlen(queryStr) - 1] = '\0';
    refStr[strlen(refStr) - 1] = '\0';
    fprintf(fp, "%d-%d\t%c\t%s\t%s\t%s\n", anno.cltTid, anno.cltIdx, strand, cltClass, queryStr, refStr);
}

/// @brief Get cluster strand
char getCltStrand(int strandFlag)
{
    char strand = '.';
    if (strandFlag > 0)
        strand = '+';
    if (strandFlag < 0)
        strand = '-';

    return strand;
}

/// @brief Generate a string which represents cluster TE-class
char *getCltClass(int numTe, int *teTable, faidx_t *teFa)
{
    static char cltClass[500] = {'\0'};
    memset(cltClass, '\0', 500);
    for (int tid = 0; tid < numTe; tid++)
    {
        if (teTable[tid] == 0)
            continue;
        strcat(cltClass, faidx_iseq(teFa, tid));
        cltClass[strlen(cltClass)] = ',';
    }

    cltClass[strlen(cltClass) - 1] = '\0';
    return cltClass;
}