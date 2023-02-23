cdef dict CLUSTER_DICT = {}
#
# ---------------------------------------------------------------
#
cdef class Cluster:
    def __init__(self):
        self.segments = []

    cpdef add(self, InsertSegment iseg):
        '''add adjacent segment into cluster'''
        self.segments.append(iseg)
#
# ---------------------------------------------------------------
#
cpdef void merge_segments(list    segments,
                          int     tid,
                          int     maxdist):
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
    
    CLUSTER_DICT[tid] = clusters
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str     fpath   = 'ex2.bam',
                         str     outpath = 'all_supp_reads.bam',
                         int32_t threads = 5,
                         uint8_t stid    = 0,
                         uint8_t maxtid  = 15,
                         uint8_t minl    = 50,
                         uint8_t maxdist = 50):
    '''build cluster
    Parameters:
    -----------
        fpath: str
            path of input BAM file.
        
        outpath: str
            file name/path for the output alignments that contain
            insert segment(s).
        
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
    cdef str  rmode = "rb"
    cdef str  wmode = "wb"
    cdef int  tid
    cdef list segments
    cdef dict SEG_DICT

    # parse alignments
    rbf = BamFile(fpath, threads, rmode)
    wbf = BamFile(outpath, threads, wmode, rbf)
    SEG_DICT = rbf.fetch(wbf, stid, maxtid, minl)

    # close file after I/O
    rbf.close(); del rbf
    wbf.close(); del wbf

    # build cluster fo each chromosome
    for tid in range(stid, maxtid):
        # sort segments by rpos
        segments = SEG_DICT[tid]
        segments.sort()

        # merge segments into cluster
        merge_segments(segments, tid, maxdist)
    
    del SEG_DICT
    return CLUSTER_DICT