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


float_t clt_entropy(int16_t nL, int16_t nM, int16_t nR, int16_t nseg) {
    float_t entropy = 0;
    if (nL > 0) entropy -= ((float_t)nL/nseg) * log2((float_t)nL/nseg);
    if (nM > 0) entropy -= ((float_t)nM/nseg) * log2((float_t)nM/nseg);
    if (nR > 0) entropy -= ((float_t)nR/nseg) * log2((float_t)nR/nseg);
    return entropy;
}


void clt_dffloc(cluster_dtype_struct clts[], int16_t nseg) {
    if (clts[0].aln1_frac > 0) clts[0].aln1_frac = clts[0].aln1_frac / nseg;
    if (clts[0].aln2_frac > 0) clts[0].aln2_frac = clts[0].aln2_frac / nseg;
    if (clts[0].aln4_frac > 0) clts[0].aln4_frac = clts[0].aln4_frac / nseg;
    if (clts[0].aln8_frac > 0) clts[0].aln8_frac = clts[0].aln8_frac / nseg;
    if (clts[0].aln16_frac > 0) clts[0].aln16_frac = clts[0].aln16_frac / nseg;
}


void cclt_feat(cluster_dtype_struct clts[], seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail, float_t div, float_t coverage, int minovh, int tid, htsFile *htsfp, bam1_t *b1, bam1_t *b2) {
    // compute cluster location flag
    clt_loc_flag(rep_ail, gap_ail, clts);

    // initilization
    uint8_t ntype=0;
    int16_t nL=0, nM=0, nR=0, nseg=0;
    for (int32_t j = clts[0].st_idx; j < clts[0].ed_idx; j++)
    {
        // ignore segment with short overhang
        if (segs[j].overhang < minovh) continue;
        
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
        nseg += 1;

        // features from TE alignment
        if (segs[j].nmap > 0)
        {
            clts[0].avg_AS += segs[j].sumAS / segs[j].nmap;
            clts[0].avg_qfrac += (float_t)segs[j].lmap / (segs[j].qed - segs[j].qst);
            clts[0].avg_div += segs[j].sumdiv / segs[j].nmap;
            // strand
            if ((segs[j].cnst & 255) > (segs[j].cnst >> 8)) {
                clts[0].strand += 1;
            } else if ((segs[j].cnst & 255) < (segs[j].cnst >> 8)) {
                clts[0].strand += 256;
            }
        } else
        {
            clts[0].avg_div += div;
        }
    }

    // single support read
    int32_t retval, second = 0;
    switch (nseg)
    {
    case 1:
        clts[0].single = 1;
        break;

    case 2:
        for (int32_t j = clts[0].st_idx; j < clts[0].ed_idx; j++)
        {
            if (segs[j].overhang < minovh) continue;
            if (second) {
                retval = bgzf_seek(htsfp->fp.bgzf, segs[j].offset, SEEK_SET);
                retval = bam_read1(htsfp->fp.bgzf, b2);
            } else {
                retval = bgzf_seek(htsfp->fp.bgzf, segs[j].offset, SEEK_SET);
                retval = bam_read1(htsfp->fp.bgzf, b1);
                second = 1;
            }
        }
        // both segments have the same qname
        retval = strcmp(bam_get_qname(b1), bam_get_qname(b2));
        if (retval == 0) clts[0].single = 1;
        break;
        
    default:
        break;
    }

    // features computation
    clts[0].nseg        = nseg / coverage;
    clts[0].ntype       = (ntype&1) + ((ntype&2)>>1) + ((ntype&4)>>2);
    clts[0].entropy     = clt_entropy(nL, nM, nR, nseg);
    clts[0].bratio      = (MIN(nL, nR) + 0.01) / (MAX(nL, nR) + 0.01);
    clts[0].lmq_frac    = clts[0].lmq_frac / nseg;
    clts[0].dclip_frac  = clts[0].dclip_frac / nseg;
    clts[0].avg_mapq    = clts[0].avg_mapq / nseg;
    clt_dffloc(clts, nseg);

    // features from TE alignment
    clts[0].avg_AS      = clts[0].avg_AS / nseg;
    clts[0].avg_qfrac   = clts[0].avg_qfrac / nseg;
    clts[0].avg_div     = (clts[0].avg_div / nseg) / div;
    // strand
    if ((clts[0].strand & 255) > (clts[0].strand >> 8)) {
        clts[0].strand = 1;
    } else if ((clts[0].strand & 255) < (clts[0].strand >> 8)) {
        clts[0].strand = 2;
    }
}