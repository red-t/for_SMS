#include "cluster_utils.h"

void clt_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, cluster_dtype_struct clts[])
{
    int32_t n = 0, d = 0x7fffffff;

    // intersecting with repeats
    query_dist_c(rep_ail, clts[0].st, clts[0].ed, 50, &n, &d);

    // intersecting with gaps
    query_dist_c(gap_ail, clts[0].st, clts[0].ed, 50, &n, &d);

    // flag of the cluster
    if (n>0){
        if (d<50) { // at boundary
            clts[0].cloc_flag = 2;
        } else { // inside repeats/gap
            clts[0].cloc_flag = 4;
        }
    } else { // at normal
        clts[0].cloc_flag = 1;
    }
}


float clt_entropy(int16_t nL, int16_t nM, int16_t nR, int16_t nseg) {
    float entropy = 0;
    if (nL > 0) entropy -= ((float)nL/nseg) * log2((float)nL/nseg);
    if (nM > 0) entropy -= ((float)nM/nseg) * log2((float)nM/nseg);
    if (nR > 0) entropy -= ((float)nR/nseg) * log2((float)nR/nseg);
    return entropy;
}


void clt_dffloc(cluster_dtype_struct clts[]) {
    if (clts[0].aln1_frac > 0) clts[0].aln1_frac = clts[0].aln1_frac / clts[0].nseg;
    if (clts[0].aln2_frac > 0) clts[0].aln2_frac = clts[0].aln2_frac / clts[0].nseg;
    if (clts[0].aln4_frac > 0) clts[0].aln4_frac = clts[0].aln4_frac / clts[0].nseg;
    if (clts[0].aln8_frac > 0) clts[0].aln8_frac = clts[0].aln8_frac / clts[0].nseg;
    if (clts[0].aln16_frac > 0) clts[0].aln16_frac = clts[0].aln16_frac / clts[0].nseg;
}


void clt_feat(cluster_dtype_struct clts[], seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail) {
    // compute cluster location flag
    clt_loc_flag(rep_ail, gap_ail, clts);

    // initilization
    uint8_t ntype=0;
    int16_t nL=0, nM=0, nR=0;
    for (int32_t j = clts[0].st_idx; j < clts[0].ed_idx; j++)
    {
        // different type of segments
        switch (segs[j].sflag)
        {
        case LEFT_CLIP:
            nL += 1;
            ntype |= LEFT_CLIP;
            break;
        case RIGHT_CLIP:
            nR += 1;
            ntype |= RIGHT_CLIP;
            break;
        case MID_INSERT:
            nM += 1;
            ntype |= MID_INSERT;
            break;

        default:
            break;
        }
        
        // short overhang
        if (segs[j].overhang < 100) clts[0].sovh_frac += 1;
        
        // low mapq
        if (segs[j].mapq < 5) clts[0].lmq_frac += 1;

        // dual-clip
        if ((segs[j].rflag & DUAL_CLIP)==5) clts[0].dclip_frac += 1;

        // segments with different loc_flag
        switch (segs[j].loc_flag)
        {
        case 1:
            clts[0].aln1_frac += 1;
            break;
        case 2:
            clts[0].aln2_frac += 1;
            break;
        case 4:
            clts[0].aln4_frac += 1;
            break;
        case 8:
            clts[0].aln8_frac += 1;
            break;
        case 16:
            clts[0].aln16_frac += 1;
            break;
        
        default:
            break;
        }
        
        // sum of mapq
        clts[0].avg_mapq += segs[j].mapq;
    }

    // features computation
    clts[0].ntype       = (ntype&1) + ((ntype&2)>>1) + ((ntype&4)>>2);
    clts[0].entropy     = clt_entropy(nL, nM, nR, clts[0].nseg);
    clts[0].bratio      = (MIN(nL, nR) + 0.01) / (MAX(nL, nR) + 0.01);
    clts[0].sovh_frac   = clts[0].sovh_frac / clts[0].nseg;
    clts[0].lmq_frac    = clts[0].lmq_frac / clts[0].nseg;
    clts[0].dclip_frac  = clts[0].dclip_frac / clts[0].nseg;
    clts[0].avg_mapq    = clts[0].avg_mapq / clts[0].nseg;
    clt_dffloc(clts);
}