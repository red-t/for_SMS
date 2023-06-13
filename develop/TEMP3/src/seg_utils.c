#include "seg_utils.h"


/***********************
 *** Segment records ***
 ***********************/

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
        if (d1<50) { // at boundary
            flag1 = 2;
        } else { // inside repeats/gap
            flag1 = 4;
        }
    } else { // at normal
        flag1 = 1;
    }

    // flag of refed
    if (n2>0){
        if (d2<50) { // at boundary
            flag2 = 2;
        } else { // inside repeats/gap
            flag2 = 4;
        }
    } else { // at normal
        flag2 = 1;
    }

    // flag of the alignment
    switch (flag1 | flag2)
    {
    // at least one side at normal region
    case 1:
    case 5:
        segs[0].loc_flag = 1;
        break;
    // both sides at repeat/gap boundary
    case 2:
        segs[0].loc_flag = 4;
        break;
    // one side at repeat/gap boundary, the other at normal region
    case 3:
        segs[0].loc_flag = 8;
        break;
    // both sides inside repeat/gap
    case 4:
        segs[0].loc_flag = 2;
        break;
    // one side at repeat/gap boundary, the other inside repeat/gap
    case 6:
        segs[0].loc_flag = 16;
        break;
    
    default:
        break;
    }
}


void cseg_feat(seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail) {
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



/*************************
 *** Alignment records ***
 *************************/

int sam_realloc_bam_data1(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    if ((bam_get_mempolicy1(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < (int)b->m_data ? b->l_data : (int)b->m_data);
            bam_set_mempolicy1(b, bam_get_mempolicy1(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}


int bam_trim1(bam1_t *src, bam1_t *dest, int32_t idx, int32_t qst, int32_t qed)
{
    // use idx as qname and compute l_qname from it
    int32_t l_seq = qed - qst;
    char qname[12];
    sprintf(qname, "%d", idx);
    int32_t l_qname = strlen(qname);

    // note: the qname is stored nul terminated and padded as described in the
    // documentation for the bam1_t struct.
    int32_t qname_nuls = 4 - l_qname % 4;

    // re-allocate the data buffer as needed.
    int32_t data_len = l_qname + qname_nuls + (l_seq + 1)/2 + 1;
    if (realloc_bam_data1(dest, data_len) < 0) {
        return -1;
    }

    dest->l_data = (int)data_len;
    dest->core.l_extranul = (uint8_t)(qname_nuls - 1);
    dest->core.flag = (src->core.flag & BAM_FREVERSE) | BAM_FUNMAP;
    dest->core.l_qname = (uint16_t)(l_qname + qname_nuls);

    // set qname
    uint8_t *cp = dest->data;
    memcpy(cp, qname, l_qname);
    int i;
    for (i = 0; i < qname_nuls; i++) {
        cp[l_qname + i] = '\0';
    }
    cp += l_qname + qname_nuls;

    // set sequence
    // 1. when qst is even, the first base of seq[qst>>1] is just the start of the query,
    // when qst is odd, the second base of seq[qst>>1] is the start of the query.
    // 2. when l_seq is even, memcpy(cp, seq, (l_seq+1)>>1) will copy sequence with length
    // of l_seq, when l_seq is odd, will copy sequence with length of (l_seq+1).
    uint8_t *seq = bam_get_seq(src);
    seq += (qst >> 1);
    if (qst & 1) { // odd qst
        if (l_seq & 1) { // odd l_seq
            memcpy(cp, seq, (l_seq+1) >> 1);
            dest->core.l_qseq = (int32_t)l_seq + 1;
        } else { // even l_seq
            memcpy(cp, seq, (l_seq+2) >> 1);
            dest->core.l_qseq = (int32_t)l_seq + 1;
        }
    } else { // even qst
        memcpy(cp, seq, (l_seq+1) >> 1);
        dest->core.l_qseq = (int32_t)l_seq;
    }

    return (int)data_len;
}