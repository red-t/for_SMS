# filename: SegmentParser.pyx

cdef inline char* getSequenceInRange(bam1_t *src, uint32_t start, uint32_t end):
    '''return python string of the sequence in a bam1_t object.'''
    cdef uint8_t *p
    cdef uint32_t k
    cdef char *s
    cdef bytes seq
    
    seq = PyBytes_FromStringAndSize(NULL, end - start)
    s   = <char*>seq
    p   = pysam_bam_get_seq(src)

    for k from start <= k < end:
        # equivalent to seq_nt16_str[bam1_seqi(s, i)] (see bam.c)
        # note: do not use string literal as it will be a python string
        s[k-start] = seq_nt16_str[p[k//2] >> 4 * (1 - k%2) & 0xf]

    return seq


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
    
    property qname:
        '''the query template name (None if not present)'''
        def __get__(self):
            cdef bam1_t * src = self._delegate
            if src.core.l_qname == 0:
                return None
            return charptr_to_str(<char *>pysam_bam_get_qname(src))
    
    property flag:
        '''properties flag'''
        def __get__(self):
            return self._delegate.core.flag
    
    property tid:
        '''tid, the index of the reference sequence'''
        def __get__(self):
            return self._delegate.core.tid
    
    property ref_start:
        '''0-based leftmost coordinate of the alignment'''
        def __get__(self):
            return self._delegate.core.pos
    
    property mapq:
        '''mapping quality'''
        def __get__(self):
            cdef bam1_t * src = self._delegate
            return pysam_get_qual(src)
    
    property qlen:
        '''
        the length of the query/read, corresponds to the length of the
        sequence supplied in the BAM/SAM file. The length of a query 
        is 0 if there is no sequence in the BAM/SAM file.
        '''
        def __get__(self):
            return self._delegate.core.l_qseq
    
    property qseq:
        '''read sequence bases'''
        def __get__(self):
            cdef bam1_t * src
            cdef str s
            src = self._delegate
            if src.core.l_qseq == 0:
                return None

            s = charptr_to_str(getSequenceInRange(src, 0, src.core.l_qseq))
            return s
    
    property is_mapped:
        '''true if read itself is mapped'''
        def __get__(self):
            return (self.flag & BAM_FUNMAP) == 0
    
    property is_forward:
        '''true if read is mapped to forward strand'''
        def __get__(self):
            return (self.flag & BAM_FREVERSE) == 0
    
    property is_secondary:
        '''true if not primary alignment'''
        def __get__(self):
            return (self.flag & BAM_FSECONDARY) != 0
    
    property is_supplementary:
        '''true if this is a supplementary alignment'''
        def __get__(self):
            return (self.flag & BAM_FSUPPLEMENTARY) != 0
    
    property reference_end:
        '''
        aligned reference position of the read on the reference genome.
        reference_end points to one past the last aligned residue.
        '''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src)
    
    property reference_length:
        '''
        aligned length of the read on the reference genome.
        equal to `reference_end - reference_start`. 
        '''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src) - self._delegate.core.pos
    
    cpdef get_tag_i(self, str tag):
        pass
    
    cpdef get_tag_f(self, str tag):
        pass
    



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
cdef int parse_cigar(bam1_t *src, list tmp_segl, uint8_t minl=50):
    '''parse CIGAR

    traverse through alignment's CIGAR, extract clip/insert segment.

    Parameters:
    -----------
        src: bam1_t
            single record pointer read from BAM file.
        tmp_segl: list
            temporary list used to stored the new InsertSegment ins
            tances
        minl: uint8_t
            minimum segment length, S/I CIGAR option with length >=
            minl will be extracted.
    '''
    cdef uint8_t orient=0
    cdef int32_t  qpos=0, rpos, qst, qend
    cdef uint32_t op, l, k, n, m
    cdef uint32_t *cigar_p
    cdef InsertSegment iseg

    n = src.core.n_cigar
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
            if l >= minl:
                qst = qpos
                qend = qpos+l
                if k>0:
                    orient = 1
                if k==m:
                    orient = 2
                iseg = makeInsertSegment(src,qst,qend,rpos,orient)
                tmp_segl.append(iseg)
            qpos += l
        elif op==2 or op==3:
            rpos += l
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