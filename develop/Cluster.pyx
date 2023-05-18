import numpy as np


CLUSTER_DTYPE = np.dtype([
    ('st',          np.int32),
    ('st_idx',      np.int32),
    ('ed',          np.int32),
    ('ed_idx',      np.int32),
    ('nseg',        np.int16),
    ('strand',      np.uint8),
    ('cloc_flag',   np.uint8),
    ('n1',          np.uint8),
    ('l1',          np.int32),
    ('n2',          np.uint8),
    ('l2',          np.int32),
    ('ntype',       np.uint8),
    ('entropy',     np.float32),
    ('bratio',      np.float32),
    ('nL',          np.int16),
    ('nM',          np.int16),
    ('nR',          np.int16),
    ('lmq_freq',    np.float32),
    ('avg_mapq',    np.float32),
])
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
    cdef int j, M
    cdef int i   = 0
    cdef int idx = 0
    cdef object clts, template
    cdef cluster_dtype_struct[::1] clts_view

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
cpdef dict build_cluster(str fpath,
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
    cdef BamFile rbf, wbf
    cdef str    rmode   = "rb"
    cdef str    wmode   = "wb"
    cdef str    outpath = "tmp.all_supp_reads.{}.bam".format(tid)
    cdef object segs, clts
    cdef dict   CLUSTER_DICT = {}
    # cdef cluster_dtype_struct[::1] clts_view

    # parse alignments & extract segments
    rbf = BamFile(fpath, threads, rmode)
    wbf = BamFile(outpath, threads, wmode, rbf)
    segs = rbf.extract_seg(wbf, tid, minl)

    # close file after I/O
    rbf.close(); del rbf
    wbf.close(); del wbf
    
    # merge segments into cluster
    segs.sort(order='rpos')
    clts = merge_segments(segs, maxdist)

    CLUSTER_DICT[tid] = (clts, segs)
    return CLUSTER_DICT