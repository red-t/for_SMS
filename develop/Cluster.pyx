import numpy as np
cdef dict CLUSTER_DICT = {}
#
# ---------------------------------------------------------------
#
cdef class Cluster:
    def __init__(self):
        self.segments = []

    cpdef add(self, InsertSegment iseg):
        '''add adjacent segment into cluster'''
        self.segments.append((iseg.qstart,
                              iseg.qend,
                              iseg.rpos,
                              iseg.stype,
                              iseg.qname))
#
# ---------------------------------------------------------------
#
cdef void merge_segments(list    segments,
                         int     tid,
                         int     maxdist,
                         object segs,
                         object clts):
    '''merge overlapped segments into cluster
    
    Parameters:
    -----------
        segments: list
            list of segments, sorted by rpos
        
        tid: int
            chromosome id of the segments
        
        maxdist: int
            max merging distance. segments with distance larger t-
            han maxdist will not be merged in to the same cluster
    '''
    cdef Cluster c
    cdef InsertSegment s
    cdef list clusters = []
    cdef int i, j, start, end, n
    
    i = 0
    n = len(segments)
    while i < n:
        s     = segments[i]
        c     = Cluster()
        start = s.rpos - 1
        end   = s.rpos + maxdist
        c.add(s)

        # try to merge segments within maxdist iteratively
        j = i + 1
        while j < n:
            s = segments[j]
            if s.rpos <= end:
                c.add(s)
                end = s.rpos + maxdist
                j += 1
            else:
                break
        
        c.start = start
        c.end   = end - maxdist
        c.nseg  = j - i
        clusters.append(c)
        # jump merged segments
        i = j
    
    CLUSTER_DICT[tid] = (clusters, segs, clts)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str     fpath,
                         int32_t threads,
                         uint8_t tid,
                         uint8_t minl,
                         uint8_t maxdist):
    '''build cluster
    Parameters:
    -----------
        fpath: str
            path of input BAM file.
        
        threads: int32_t
            nummber of threads used for BAM file I/O.
        
        stid: uint8_t
            start tid, cluster building start from stid-th contig.
        
        maxtid: uint8_t
            max tid, cluster building end at (maxtid-1)-th contig.
        
        minl: uint8_t
            minimun segment length, cigar operation with length < 
            minl will not be used to create insert segment.
        
        maxdist: uint8_t
            max merging distance, segments with distance larger t-
            han maxdist will not be merged in to the same cluster.
    
    Returns:
    --------
        CLUSTER_DICT: dict
            dictionary of clusters, tid -> list_of_Cluster.
    '''
    cdef BamFile rbf, wbf
    cdef list segments
    cdef dict SEG_DICT
    cdef str  rmode   = "rb"
    cdef str  wmode   = "wb"
    cdef str  outpath = "tmp.all_supp_reads.{}.bam".format(tid)
    cdef object clts  = 0

    # parse alignments
    rbf = BamFile(fpath, threads, rmode)
    wbf = BamFile(outpath, threads, wmode, rbf)
    SEG_DICT = rbf.extract_seg(wbf, tid, minl)

    # close file after I/O
    rbf.close(); del rbf
    wbf.close(); del wbf
    
    # merge segments into cluster
    segments, segs = SEG_DICT[tid]
    segments.sort()
    merge_segments(segments, tid, maxdist, segs, clts)
    
    del SEG_DICT
    return CLUSTER_DICT





#
# ---------------------------------------------------------------
#
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
cdef object merge_segments1(seg_dtype_struct[::1] segs,
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
    M = clts_view.shape[0] - 10 
    
    while i < segs.shape[0]:
        if idx == M:
            clts = np.concatenate((clts, template))
            clts_view = clts
            M = clts_view.shape[0] - 10 

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
cpdef dict build_cluster1(str fpath,
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
        CLUSTER_DICT: dict
            dictionary of clusters, tid -> list_of_Cluster.
    '''
    cdef BamFile rbf, wbf
    cdef list   segments
    cdef dict   SEG_DICT
    cdef str    rmode   = "rb"
    cdef str    wmode   = "wb"
    cdef str    outpath = "tmp.all_supp_reads.{}.bam".format(tid)
    cdef object segs, clts
    # cdef cluster_dtype_struct[::1] clts_view

    # parse alignments & extract segments
    rbf = BamFile(fpath, threads, rmode)
    wbf = BamFile(outpath, threads, wmode, rbf)
    SEG_DICT = rbf.extract_seg(wbf, tid, minl)

    # close file after I/O
    rbf.close(); del rbf
    wbf.close(); del wbf
    
    # merge segments into cluster
    segments, segs = SEG_DICT[tid]
    # new version
    segs.sort(order='rpos')
    clts = merge_segments1(segs, maxdist)
    # for comparison
    segments.sort()
    merge_segments(segments, tid, maxdist, segs, clts)
    
    del SEG_DICT
    return CLUSTER_DICT