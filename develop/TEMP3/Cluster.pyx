import numpy as np
from subprocess import Popen, DEVNULL


SEG_DTYPE = np.dtype([
    ('flag',        np.uint16),
    ('mapq',        np.uint8),
    ('qst',         np.int32),
    ('qed',         np.int32),
    ('rpos',        np.int32),
    ('lqseq',       np.int32),
    ('sflag',       np.uint8),
    ('rflag',       np.uint8),
    ('offset',      np.int64),
    ('refst',       np.int32),
    ('refed',       np.int32),
    ('ith',         np.uint8),
    ('nseg',        np.uint8),
    ('overhang',    np.int32),
    ('nmatch',      np.int32),
    ('loc_flag',    np.uint8),
])


TEALN_DTYPE = np.dtype([
    ('idx',     np.int32),
    ('AS',      np.int32),
    ('qst',     np.int32),
    ('qed',     np.int32),
    ('flag',    np.int16),
])


CLUSTER_DTYPE = np.dtype([
    ('st',          np.int32),
    ('ed',          np.int32),
    ('st_idx',      np.int32),
    ('ed_idx',      np.int32),
    ('nseg',        np.int16),
    ('strand',      np.uint8),
    ('cloc_flag',   np.uint8),
    ('ntype',       np.uint8),
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
cdef object extract_seg(Iterator ite, int minl):
    cdef:
        int32_t     retval, M, N=0
        seg_dtype_struct[::1]   segs_view

    segs = np.zeros(10000, dtype=SEG_DTYPE)
    template  = np.zeros(10000, dtype=SEG_DTYPE)
    segs_view = segs
    M = segs_view.shape[0] - 20
    
    while 1:
        retval = ite.cnext1()
        if retval > 0:
            if ite.b.core.n_cigar == 0:
                continue
            if N > M:
                segs = np.concatenate((segs, template))
                segs_view = segs
                M = segs_view.shape[0] - 20

            retval = parse_cigar(ite.b, &segs_view[N], ite.offset, minl)
            N += retval
            continue

        del template
        return segs[:N]
#
# ---------------------------------------------------------------
#
cdef seg_feat(seg_dtype_struct[::1] segs,
              ailist_t *rep_ail,
              ailist_t *gap_ail):
    cdef ssize_t i

    for i in range(segs.shape[0]):
        cseg_feat(&segs[i], rep_ail, gap_ail)
#
# ---------------------------------------------------------------
#
cdef trim_seg(int tid,
              int threads,
              BamFile rbf,
              Iterator ite,
              seg_dtype_struct[::1] segs):
    cdef:
        str outpath = "tmp.all_supp_reads.{}.fa".format(tid)
        BamFile wbf = BamFile(outpath, threads, "wF", rbf)
        bam1_t *dest = bam_init1()
        int i, retval
    
    for i in range(segs.shape[0]):
        if aln_is_second(segs[i].flag):
            continue
        retval = ite.cnext3(segs[i].offset)
        if retval < 0:
            raise StopIteration
        # trim alignment by bam_trim1, and write out
        bam_trim1(ite.b, dest, i, segs[i].qst, segs[i].qed)
        wbf.write(dest)
    
    bam_destroy1(dest); wbf.close(); del wbf
#
# ---------------------------------------------------------------
#
cdef align_mm2(int tid,
               int threads,
               str ref,
               str preset):
    cdef:
        int retval
        str cmd_mm2 = "minimap2 -t {} -aYx {} {} tmp.all_supp_reads.{}.fa | " \
                      "samtools view -@ {} -bhS -o tmp.all_supp_reads.{}.bam -".format(threads, preset, ref, tid, threads, tid)

    proc = Popen([cmd_mm2], stderr=DEVNULL, shell=True, executable='/bin/bash')
    retval = proc.wait()
    if retval != 0:
        raise Exception("Error: minimap2 failed for tid: {}".format(tid))
#
# ---------------------------------------------------------------
#
cdef object extract_tealn(Iterator ite):
    cdef:
        int32_t     retval, M, N=0
        tealn_dtype_struct[::1]   tealns_view

    tealns = np.zeros(10000, dtype=TEALN_DTYPE)
    template  = np.zeros(10000, dtype=TEALN_DTYPE)
    tealns_view = tealns
    M = tealns_view.shape[0] - 20
    
    while 1:
        retval = ite.cnext2()
        if retval > 0:
            if bam_filtered(ite.b):
                continue
            if N > M:
                tealns = np.concatenate((tealns, template))
                tealns_view = tealns
                M = tealns_view.shape[0] - 20

            parse_tealns(ite.b, &tealns_view[N])
            N += 1
            continue

        del template
        return tealns[:N]
#
# ---------------------------------------------------------------
#
cdef object seg_feat_te(seg_dtype_struct[::1] segs, int tid, int threads):
    cdef:
        BamFile  rbf = BamFile("tmp.all_supp_reads.{}.bam".format(tid), threads, "rb")
        Iterator ite = Iterator(rbf)
        object   alns
        tealn_dtype_struct[::1] tealns_view
    
    alns = extract_tealn(ite)
    alns.sort(order=['idx', 'qst'])
    tealns_view = alns

    del ite; rbf.close(); del rbf
    return alns
#
# ---------------------------------------------------------------
#
cdef object merge_seg(seg_dtype_struct[::1] segs,
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
cdef clt_feat(cluster_dtype_struct[::1] clts,
              seg_dtype_struct[::1] segs,
              ailist_t *rep_ail,
              ailist_t *gap_ail):
    cdef:
        ssize_t i, j
    
    for i in range(clts.shape[0]):
        cclt_feat(&clts[i], &segs[0], rep_ail, gap_ail)
    
    # TO DO:
    # 1. refine cclt_feat to compute clt features related to TE alignments
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str fpath,
                         str rep_path,
                         str gap_path,
                         str teref,
                         str preset,
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
    ##############################################
    ### 1. parse alignments & extract segments ###
    ##############################################
    cdef:
        object   segs
        BamFile  rbf = BamFile(fpath, threads, "rb")
        Iterator ite = Iterator(rbf, tid)
    
    # extract segments
    segs = extract_seg(ite, minl)

    ############################################
    ### 2. compute features for each segment ###
    ############################################
    cdef:
        bytes repfn = rep_path.encode()
        bytes gapfn = gap_path.encode()
        const char *chrom = sam_hdr_tid2name(rbf.hdr, tid)
        ailist_t *rep_ail = ailist_init()
        ailist_t *gap_ail = ailist_init()
        seg_dtype_struct[::1] segs_view = segs

    # construct AIList
    readBED(rep_ail, repfn, chrom); ailist_construct(rep_ail, 20)
    readBED(gap_ail, gapfn, chrom); ailist_construct(gap_ail, 20)

    # compute features for each segment
    seg_feat(segs_view, rep_ail, gap_ail)

    # trimmed & write out segment sequence
    segs.sort(order='rpos')
    trim_seg(tid, threads, rbf, ite, segs_view)
    rbf.close(); del rbf; del ite

    # TO DO:
    # 1. write a function that map the segment sqeunces to TE CSS (by calling minimap2)
    # 2. write a function to read TE alignments & compute new features for segments
    # align segment sequences to TE CSS
    align_mm2(tid, threads, teref, preset)
    alns = seg_feat_te(segs_view, tid, threads)
    
    ############################
    ### 3. construct cluster ###
    ############################
    cdef:
        object clts
        dict   CLUSTER_DICT = {}

    # merge segments into cluster
    clts = merge_seg(segs, maxdist)

    ############################################
    ### 4. compute features for each cluster ###
    ############################################
    cdef cluster_dtype_struct[::1] clts_view = clts
    clt_feat(clts_view, segs_view, rep_ail, gap_ail)

    # free AIList
    ailist_destroy(rep_ail); ailist_destroy(gap_ail)

    CLUSTER_DICT[tid] = (clts, segs, alns)
    return CLUSTER_DICT