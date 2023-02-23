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
cpdef void merge_segments(BamFile wbf,
                          list    segments,
                          int     tid,
                          int     fsize):
    '''merge overlapped segments into cluster'''
    cdef list clusters = []
    cdef Cluster c
    cdef InsertSegment s
    cdef int i, j, start, end, n
    
    # 开始合并
    i = 0
    n = len(segments)
    while i < n:
        s     = segments[i]
        c     = Cluster()
        start = s.rpos - 1
        end   = s.rpos + fsize

        # 添加 cluster 的第一个 segment
        c.add(s)
        wbf.write(s)

        # 寻找需要合并的区间
        j = i + 1
        while j < n:
            s = segments[j]
            if s.rpos <= end:
                # 如果有 overlap，将segment添加到cluster当中
                c.add(s)
                wbf.write(s)
                # 更新区间的边界
                end = s.rpos + fsize
                j += 1
            else:
                break
        
        # 将 cluster 添加到结果当中
        c.start = start
        c.end   = end - fsize
        c.nseg  = j - i
        clusters.append(c)
        # 跳过已经合并的区间
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
                         uint8_t fsize   = 50):
    '''build cluster'''
    cdef BamFile rbf, wbf
    cdef str  rmode = "rb"
    cdef str  wmode = "wb"
    cdef int  tid
    cdef list segments
    cdef dict SEG_DICT

    # read and parse alignments
    rbf = BamFile(fpath, threads, rmode)
    SEG_DICT = rbf.fetch(stid, maxtid, minl)

    # create a BamFile object for writing
    wbf = BamFile(outpath, threads, wmode, rbf)
    rbf.close(); del rbf

    # call merge_intervals() for segments on each chromosome 
    for tid in range(stid, maxtid):
        # 按照 rpos 排序
        segments = SEG_DICT[tid]
        segments.sort()

        # 合并相邻的区间，创建 cluster
        merge_segments(wbf, segments, tid, fsize)
    
    wbf.close(); del wbf; del SEG_DICT
    return CLUSTER_DICT