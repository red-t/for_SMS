#include "io_utils.h"

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

int outputRefFlankSeqs(char *refFn, Cluster *cltArray, int startIdx, int endIdx)
{
    hts_pos_t seqLen;
    int leftStart, leftEnd, rightStart, rightEnd;
    faidx_t *refFa = fai_load((const char*)refFn);
    char *outFn = (char *)malloc(50 * sizeof(char));

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

    free(outFn);
    fai_destroy(refFa);
    return 0;
}