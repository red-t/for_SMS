#include "seg_utils.h"

void aln_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, seg_dtype_struct segs[])
{
    int32_t n1 = 0, n2 = 0;
    int32_t d1 = 0x7fffffff, d2 = 0x7fffffff;
    uint8_t flag1, flag2;

    // intersecting with repeats
    query_dist_p(rep_ail, segs[0].refst, 50, &n1, &d1);
    query_dist_p(rep_ail, segs[0].refed, 50, &n2, &d2);

    // intersecting with gaps
    query_dist_p(gap_ail, segs[0].refst, 50, &n1, &d1);
    query_dist_p(gap_ail, segs[0].refed, 50, &n2, &d2);

    // flag of refst
    if (n1>0){
        if (d1<50) {
            // at boundary
            flag1 = 1;
        } else {
            // inside repeats/gap
            flag1 = 2;
        }
    } else {
        // at normal
        flag1 = 0;
    }

    // flag of refed
    if (n2>0){
        if (d2<50) {
            // at boundary
            flag2 = 1;
        } else {
            // inside repeats/gap
            flag2 = 2;
        }
    } else {
        // at normal
        flag2 = 0;
    }

    // flag of the alignment
    if (flag1==1 && flag2==1){ // both sides at repeat/gap boundary
        segs[0].loc_flag = 4;
    } else if (flag1==2 && flag2==2) { // both sides inside repeat/gap
        segs[0].loc_flag = 2;
    } else if (flag1==1 || flag2==1) { // one side at repeat/gap boundary
        if (flag1==0 || flag2==0) { // the other side at normal region
            segs[0].loc_flag = 8;
        } else { // the other side inside repeat/gap
            segs[0].loc_flag = 16;
        }
    } else { // at least one side at normal region
        segs[0].loc_flag = 1;
    }
}


void seg_feat(seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail) {
    // update overhang
    if (segs[0].sflag & LEFT_CLIP) {
        segs[0].overhang = segs[0].nmatch;
    } else if (segs[0].sflag & MID_INSERT) {
        segs[0].overhang = update_overhang(segs[0].overhang, segs[0].nmatch);
    }

    // compute alignment location flag
    if (segs[0].ith == 0) {
        aln_loc_flag(rep_ail, gap_ail, segs);
        if (segs[0].nseg > 1) {
            for (uint8_t j = 1; j < segs[0].nseg; j++)
            {
                segs[j].loc_flag = segs[0].loc_flag;
            }
        }
    }
    
    // TO DO
    // 1. compute the "trimmed" start & end for each segment (need a C function)
    // 2. read alignment with specified offset by seek (need a function)
    // 3. modified & write alignment by bam_set_qname, bam_set1, sam_write1 (need a C function)
}