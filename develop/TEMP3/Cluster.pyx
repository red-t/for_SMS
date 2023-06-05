import numpy as np


CLUSTER_DTYPE = np.dtype([
    ('st',          np.int32),
    ('ed',          np.int32),
    ('st_idx',      np.int32),
    ('ed_idx',      np.int32),
    ('nseg',        np.int16),
    ('strand',      np.uint8),
    ('cloc_flag',   np.uint8),
    ('ntype',       np.uint8),
    ('nL',          np.int16),
    ('nM',          np.int16),
    ('nR',          np.int16),
    ('entropy',     np.float32),
    ('bratio',      np.float32),
    ('sovh_frac',   np.float32),
    ('lmq_frac',    np.float32),
    ('dclip_frac',  np.float32),
    ('aln1_frac',   np.float32),
    ('aln2_frac',   np.float32),
    ('aln4_frac',   np.float32),
    ('aln8_frac',   np.float32),
    ('aln16_frac',  np.float32),
    ('avg_mapq',    np.float32),
])
#
# ---------------------------------------------------------------
#
cdef object extract_seg(Iterator ite):
    cdef ssize_t i
    for i in ite:
        continue

    return ite.segs[:ite.N,]
#
# ---------------------------------------------------------------
#
cdef compute_seg_feat(seg_dtype_struct[::1] segs,
                      ailist_t *rep_ail,
                      ailist_t *gap_ail):
    cdef ssize_t i, j

    for i in range(segs.shape[0]):
        # update overhang
        if segs[i].sflag & LEFT_CLIP:
            segs[i].overhang = segs[i].nmatch
        elif segs[i].sflag & MID_INSERT:
            segs[i].overhang = update_overhang(segs[i].overhang, segs[i].nmatch)

        # compute alignment location flag
        if segs[i].ith == 0:
            aln_loc_flag(rep_ail, gap_ail, &segs[i])
            if segs[i].nseg > 1:
                for j in range(1, segs[i].nseg):
                    segs[i+j].loc_flag = segs[i].loc_flag
        
        # TO DO
        # 1. compute the "trimmed" start & end for each segment (need a C function)
        # 2. read alignment with specified offset by seek (need a function)
        # 3. modified & write alignment by bam_set_qname, bam_set1, sam_write1 (need a C function)
#
# ---------------------------------------------------------------
#
cdef object merge_segments(seg_dtype_struct[::1] segs,
                                            int maxdist):
    '''merge overlapped segments into cluster
    
    Parameters:
    -----------
        segs: 
            typed memoryview of a structed numpy arrary, which stores
            features of the extracted segments.

        maxdist: int
            max merging distance. segments with distance larger t-
            han maxdist will not be merged in to the same cluster
    '''
    cdef:
        int j, M
        int i   = 0
        int idx = 0
        object clts, template
        cluster_dtype_struct[::1] clts_view

    clts      = np.zeros(10000, dtype=CLUSTER_DTYPE)
    template  = np.zeros(10000, dtype=CLUSTER_DTYPE)
    clts_view = clts
    M = clts_view.shape[0] - 20 
    
    while i < segs.shape[0]:
        if idx > M:
            clts = np.concatenate((clts, template))
            clts_view = clts
            M = clts_view.shape[0] - 20

        # initialize the idx-th cluster with the first segment
        clts_view[idx].st     = segs[i].rpos - 1
        clts_view[idx].st_idx = i
        clts_view[idx].ed     = segs[i].rpos + maxdist

        # try to merge segments within maxdist iteratively
        j = i + 1
        while j < segs.shape[0]:
            # next segment's rpos VS current cluster end
            if segs[j].rpos <= clts_view[idx].ed:
                clts_view[idx].ed = segs[j].rpos + maxdist
                j += 1
            else:
                break
        
        # update cluster's ed, ed_idx, nseg with the last segment
        # idx-th cluster contains segments in [i,j); nseg = j-i
        clts_view[idx].ed     = clts_view[idx].ed - maxdist
        clts_view[idx].ed_idx = j
        clts_view[idx].nseg   = j - i

        # jump merged segments
        i = j
        idx += 1
    
    del template
    return clts[:idx]
#
# ---------------------------------------------------------------
#
cdef compute_clt_feat1(cluster_dtype_struct[::1] clts,
                      seg_dtype_struct[::1] segs,
                      ailist_t *rep_ail,
                      ailist_t *gap_ail):
    cdef:
        ssize_t i, j
    
    for i in range(clts.shape[0]):
        # compute cluster location flag
        clt_loc_flag(rep_ail, gap_ail, &clts[i])

        for j in range(clts[i].st_idx, clts[i].ed_idx):
            # different type of segments
            if segs[j].sflag & LEFT_CLIP:
                clts[i].nL += 1
                clts[i].ntype |= LEFT_CLIP
            elif segs[j].sflag & RIGHT_CLIP:
                clts[i].nR += 1
                clts[i].ntype |= RIGHT_CLIP
            else:
                clts[i].nM += 1
                clts[i].ntype |= MID_INSERT
            
            # short overhang
            if segs[j].overhang < 100:
                clts[i].sovh_frac += 1
            
            # low mapq
            if segs[j].mapq < 5:
                clts[i].lmq_frac += 1
            
            # dual-clip
            if (segs[j].rflag & DUAL_CLIP)==5:
                clts[i].dclip_frac += 1
            
            # algnments with different loc_flag
            if segs[j].loc_flag & 1:
                clts[i].aln1_frac += 1
            elif segs[j].loc_flag & 2:
                clts[i].aln2_frac += 1
            elif segs[j].loc_flag & 4:
                clts[i].aln4_frac += 1
            elif segs[j].loc_flag & 8:
                clts[i].aln8_frac += 1
            elif segs[j].loc_flag & 16:
                clts[i].aln16_frac += 1
            
            # sum of mapq
            clts[i].avg_mapq += segs[j].mapq

        clt_ntype(&clts[i]) # number of segment types
        clt_entropy(&clts[i]) # compute entropy
        clt_bratio(&clts[i]) # compute bratio
        clt_sovh(&clts[i]) # compute short overhang fraction
        clt_lmq(&clts[i]) # compute low mapq fraction
        clt_dclip(&clts[i]) # compute dual-clip fraction
        clt_dffloc(&clts[i]) # compute fraction of alignment with different loc_flag
        clt_avgmapq(&clts[i]) # compute average mapqs


# cdef compute_clt_feat(cluster_dtype_struct[::1] clts,
#                       seg_dtype_struct[::1] segs,
#                       ailist_t *rep_ail,
#                       ailist_t *gap_ail):
#     cdef:
#         ssize_t i, j
    
#     for i in range(clts.shape[0]):
#         clt_feat(&clts[i], &segs[0], rep_ail, gap_ail)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str fpath,
                         str rep_path,
                         str gap_path,
                         int threads,
                         int tid,
                         int minl,
                         int maxdist):
    '''build cluster
    Parameters:
    -----------
        fpath: str
            path of input BAM file.
        
        threads: int
            nummber of threads used for BAM file I/O.
        
        tid: int
            tid of the target chromosome.
        
        minl: int
            minimun segment length, cigar operation with length < 
            minl will not be used to create insert segment.
        
        maxdist: int
            max merging distance, segments with distance larger t-
            han maxdist will not be merged in to the same cluster.
    
    Returns:
    --------
        clts: object
            strctured numpy array with `dtype=CLUSTER_DTYPE`.
    '''
    ### 1. parse alignments & extract segments ###
    cdef:
        object  segs
        BamFile rbf, wbf
        str     outpath = "tmp.all_supp_reads.{}.bam".format(tid)

    rbf  = BamFile(fpath, threads, "rb")
    wbf  = BamFile(outpath, threads, "wb", rbf)
    ite  = Iterator(rbf, wbf, tid, minl)
    segs = extract_seg(ite)

    ### 2. compute features for segments ###
    cdef:
        bytes repfn = rep_path.encode()
        bytes gapfn = gap_path.encode()
        const char *chrom = sam_hdr_tid2name(rbf.hdr, tid)
        ailist_t *rep_ail = ailist_init()
        ailist_t *gap_ail = ailist_init()
        seg_dtype_struct[::1] segs_view = segs

    # construct AIList
    readBED(rep_ail, repfn, chrom)
    readBED(gap_ail, gapfn, chrom)
    ailist_construct(rep_ail, 20)
    ailist_construct(gap_ail, 20)

    # compute features for each segment
    compute_seg_feat(segs_view, rep_ail, gap_ail)

    # close file after I/O
    rbf.close(); del rbf
    wbf.close(); del wbf
    del ite
    
    ### 3. construct cluster ###
    cdef:
        object clts
        dict   CLUSTER_DICT = {}

    # merge segments into cluster
    segs.sort(order='rpos')
    clts = merge_segments(segs, maxdist)

    # compute features for each cluster
    cdef cluster_dtype_struct[::1] clts_view = clts
    compute_clt_feat1(clts_view, segs_view, rep_ail, gap_ail)

    # free AIList
    ailist_destroy(rep_ail)
    ailist_destroy(gap_ail)

    CLUSTER_DICT[tid] = (clts, segs)
    return CLUSTER_DICT