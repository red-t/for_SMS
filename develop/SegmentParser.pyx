# filename: SegmentParser.pyx

cdef class InsertSegment:
    def __init__(self):
        # see bam_init1
        self._delegate = <bam1_t*>calloc(1, sizeof(bam1_t))
        if self._delegate == NULL:
            raise MemoryError("could not allocated memory of {} bytes".format(sizeof(bam1_t)))
        # allocate some memory. If size is 0, calloc does not return a
        # pointer that can be passed to free() so allocate 40 bytes
        # for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc(self._delegate.m_data, 1)
        if self._delegate.data == NULL:
            raise MemoryError("could not allocate memory of {} bytes".format(self._delegate.m_data))
        self._delegate.l_data = 0
        # set some data to make read approximately legit.
        # Note, SAM writing fails with q_name of length 0
        self._delegate.core.l_qname = 0
        self._delegate.core.tid = -1
        self._delegate.core.pos = -1
        self._delegate.core.mtid = -1
        self._delegate.core.mpos = -1

    def __dealloc__(self):
        bam_destroy1(self._delegate)
#
# ---------------------------------------------------------------
#
cdef InsertSegment makeInsertSegment(bam1_t *src, int32_t qstart, int32_t qend,
                                     int32_t rpos, uint8_t orient):
    '''return an InsertSegment object constructed from `src`
    
    Parameters:
    -----------
        src: bam1_t
            a bam1_t pointer, point to a single alignment record.
        qstart: int32_t
            query start of the insert/clip segment, include.
        qend: int32_t
            query end of the insert/clip segment, not include.
        rpos: int32_t
            reference position of the insert/clip segment, corresponds to `qend`.
        orient: uint8_t
            orient of the segment, 0:left-clip, 1:mid-insert, 2:right-clip
    
    Returns:
    --------
        dest: InsertSegment
            a InsertSegment object, with independent copy of `src`.
    '''
    # note that the following does not call __init__
    cdef InsertSegment dest = InsertSegment.__new__(InsertSegment)
    dest._delegate = bam_dup1(src)
    dest.qstart = qstart
    dest.qend = qend
    dest.rpos = rpos
    dest.orient = orient
    return dest
#
# ---------------------------------------------------------------
#
cdef list parse_cigar(bam1_t *src, uint8_t minl=49):
    '''parse CIGAR

    traverse through alignment's CIGAR, extract clip/insert segment.

    Parameters:
    -----------
        src: bam1_t
            single record pointer read from BAM file.
        minl: uint8_t
            minimum segment length, S/I CIGAR option with length>minl will be extracted.

    Returns:
    --------
        segl: list
            list of InsertSegment instances.
    '''
    cdef uint8_t orient=0
    cdef int32_t  qpos=0, rpos, qst, qend
    cdef uint32_t op, l, k, n, m
    cdef uint32_t *cigar_p
    cdef list segl=[]
    cdef InsertSegment iseg

    n = src.core.n_cigar
    if n == 0:
        return segl
    
    m = n-1
    cigar_p = bam_get_cigar(src)
    rpos = src.core.pos
    for k from 0 <= k < n:      # read through alignment's CIGAR
        op = cigar_p[k] & BAM_CIGAR_MASK
        l = cigar_p[k] >> BAM_CIGAR_SHIFT
        if op==0 or op==7 or op==8:
            qpos += l
            rpos += l
        elif op==1 or op==4 or op==6:
            if l > minl:
                qst = qpos
                qend = qpos+l
                if k>0:
                    orient = 1
                if k==m:
                    orient = 2
                iseg = makeInsertSegment(src,qst,qend,rpos,orient)
                segl.append(iseg)
            qpos += l
        elif op==2 or op==3:
            rpos += l

    return segl
#
# ---------------------------------------------------------------
#
cdef traverse_alignment():
    pass
#
# ---------------------------------------------------------------
#
cdef merge_intervals():
    pass
#
# ---------------------------------------------------------------
#
cdef build_cluster():
    pass