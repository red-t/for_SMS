#include "anno_utils.h"

/// @brief init single annotate record
void initAnno(bam1_t *bamRecord, Anno *anno, int idx)
{
    int numCigar = bamRecord->core.n_cigar;
    uint32_t *cigarArray = bam_get_cigar(bamRecord);

    int queryStart = 0;
    if (bam_cigar_op(cigarArray[0]) == BAM_CSOFT_CLIP)
        queryStart = bam_cigar_oplen(cigarArray[0]);

    int queryEnd = bamRecord->core.l_qseq;
    if (bam_cigar_op(cigarArray[numCigar-1]) == BAM_CSOFT_CLIP)
        queryEnd -= bam_cigar_oplen(cigarArray[numCigar-1]);

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

    if (bam_is_rev(bamRecord)) {
        anno->queryStart = bamRecord->core.l_qseq - queryEnd;
        anno->queryEnd = bamRecord->core.l_qseq - queryStart;
    }
}

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
    while (1)
    {
        int returnValue = bam_read1(inputBam->fp.bgzf, bamRecord);
        if (returnValue < 0)
            break;
        if (bamIsInvalid(bamRecord))
            continue;

        initAnno(bamRecord, &annoArray[numAnno], idx);
        leftMost = (annoArray[numAnno].queryStart < leftMost) ? annoArray[numAnno].queryStart : leftMost;
        rightMost = (annoArray[numAnno].queryEnd > rightMost) ? annoArray[numAnno].queryEnd : rightMost;
        numAnno++;
    }

    if (numAnno > 0)
        cluster->flag |= CLT_TE_MAP;·

    // annoPolyA

    if (bamRecord != NULL) {bam_destroy1(bamRecord); bamRecord=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    if (inputFn != NULL) {free(inputFn); inputFn=NULL;}

    return numAnno;
}
