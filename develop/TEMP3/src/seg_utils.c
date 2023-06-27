#include "seg_utils.h"


/***********************
 *** Segment records ***
 ***********************/

int parse_cigar(bam1_t *bam, seg_dtype_struct segs[], int64_t offset, int minl)
{   
    uint32_t *cigar = bam_get_cigar(bam);
    uint8_t  idx = 0;
    uint8_t  rflag = 0;
    int nmatch = 0;
    int qpos = 0;
    int rpos = bam->core.pos;
    int n = bam->core.n_cigar;
    int m = n-1;

    // traverse alignment's CIGAR
    for (int i = 0; i < n; i++)
    {
        int len = bam_cigar_oplen(cigar[i]);
        switch (bam_cigar_op(cigar[i]))
        {
        // match
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            qpos += len;
            rpos += len;
            nmatch += len;
            break;
        // clip or insert
        case BAM_CSOFT_CLIP:
        case BAM_CINS:
            if (len >= minl)
            {
                if (i == 0)
                {
                    segs[idx].qst   = qpos;
                    segs[idx].qed   = qpos + len;
                    segs[idx].rpos  = rpos;
                    segs[idx].sflag = LEFT_CLIP;
                    segs[idx].ith   = idx;
                    segs[idx].overhang = nmatch;
                    rflag |= LEFT_CLIP;
                } else if (i == m)
                {
                    segs[idx].qst   = qpos;
                    segs[idx].qed   = qpos + len;
                    segs[idx].rpos  = rpos;
                    segs[idx].sflag = RIGHT_CLIP;
                    segs[idx].ith   = idx;
                    segs[idx].overhang = nmatch;
                    rflag |= RIGHT_CLIP;
                } else
                {
                    segs[idx].qst   = qpos;
                    segs[idx].qed   = qpos + len;
                    segs[idx].rpos  = rpos;
                    segs[idx].sflag = MID_INSERT;
                    segs[idx].ith   = idx;
                    segs[idx].overhang = nmatch;
                    rflag |= MID_INSERT;
                }
                idx += 1;
            }
            qpos += len;
            break;
        // del or skip
        case BAM_CDEL:
        case BAM_CREF_SKIP:
            rpos += len;
            break;
        // defalut
        default:
            break;
        }
    }

    // the same attributes of all segments
    if (idx > 0)
    {
        for (uint8_t i = 0; i < idx; i++)
        {
            segs[i].flag    = bam->core.flag;
            segs[i].mapq    = bam->core.qual;
            segs[i].rflag   = rflag;
            segs[i].offset  = offset;
            segs[i].refst   = bam->core.pos;
            segs[i].refed   = rpos;
            segs[i].nseg    = idx;
            segs[i].nmatch  = nmatch;
        }
    }

    return (int)idx;
}


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


void cseg_feat(seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail)
{
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
}


void cseg_feat_te(seg_dtype_struct segs[], tealn_dtype_struct tealns[], int i)
{
    int idx = tealns[i].idx;

    // compute mapped query length
    if (segs[idx].nmap == 0)
    {
        segs[idx].lmap += tealns[i].qed - tealns[i].qst;
    } else
    { // i-th and (i-1)-th alignment come from the same segment
        int j = i - 1;
        if (tealns[i].qst < tealns[j].qed) { // i-th alignment overlap with (i-1)-th alignment
            if (tealns[i].qed > tealns[j].qed) {
                segs[idx].lmap += tealns[i].qed - tealns[j].qed;
            } else { // i-th alignment covered by (i-1)-th alignment
                tealns[i].qed = tealns[j].qed;
                return;
            }
        }
        else { // no overlap
            segs[idx].lmap += tealns[i].qed - tealns[i].qst;
        }
    }

    // compute average alignment score, nmap, consistency
    segs[idx].sumAS += (float_t)tealns[i].AS/(tealns[i].qed - tealns[i].qst);
    segs[idx].sumdiv += tealns[i].div;
    segs[idx].nmap += 1;
    if ((segs[idx].flag & BAM_FREVERSE) == (tealns[i].flag & BAM_FREVERSE)) {
        segs[idx].cnst += 1;
    } else {
        segs[idx].cnst += (1 << 8);
    }
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


void origion_qpos(int32_t *qst, int32_t *qed, bam1_t *bam)
{
    int32_t st, ed;
    uint32_t *cigar = bam_get_cigar(bam);
    int32_t n = bam->core.n_cigar;

    if (bam_is_rev(bam))
    { // bam record is reverse
        if (bam_cigar_op(cigar[n-1]) == BAM_CSOFT_CLIP)
        {
            st = bam_cigar_oplen(cigar[n-1]);
        } else
        {
            st = 0;
        }

        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
        {
            ed = bam->core.l_qseq - bam_cigar_oplen(cigar[0]);
        } else
        {
            ed = bam->core.l_qseq;
        }
    } else { // bam record is forward
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
        {
            st = bam_cigar_oplen(cigar[0]);
        } else
        {
            st = 0;
        }

        if (bam_cigar_op(cigar[n-1]) == BAM_CSOFT_CLIP)
        {
            ed = bam->core.l_qseq - bam_cigar_oplen(cigar[n-1]);
        } else
        {
            ed = bam->core.l_qseq;
        }
    }

    *qst = st;
    *qed = ed;
}


void parse_tealns(bam1_t *bam, tealn_dtype_struct tealns[])
{   int32_t qst, qed;

    origion_qpos(&qst, &qed, bam);
    tealns[0].idx = atoi(bam_get_qname(bam));
    tealns[0].AS = bam_aux2i(bam_aux_get(bam, "AS"));
    tealns[0].qst = qst;
    tealns[0].qed = qed;
    tealns[0].div = bam_aux2f(bam_aux_get(bam, "de"));
    tealns[0].flag = bam->core.flag;
}