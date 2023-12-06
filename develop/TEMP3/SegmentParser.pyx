# filename: SegmentParser.pyx

######################
### CIGAR resolver ###
######################

cdef int parse_cigar(bam1_t *src,
                     seg_dtype_struct[::1] segs,
                     int32_t N,
                     int64_t offset,
                     int32_t minl=50):
    '''parse CIGAR

    traverse through alignment's CIGAR, extract clip/insert segment.

    Parameters:
    -----------
        src: bam1_t*
            pointer to a single alignment record in BAM file.

        segs: seg_dtype_struct[::1]
            typed memoryview of a structed numpy arrary, which will 
            be used to store features of extracted segments.

        N: int32_t
            cumulative total number of extracted segments.
        
        offset: int64_t
            offset of current record in the BAM file.
        
        minl: int32_t
            minimum segment length, S/I CIGAR option with length >=
            minl will be extracted.
    
    Returns:
    --------
        nseg: int8_t
            number of insert segments extracted from this alignment
    '''
    cdef int8_t         ith         = 0
    cdef int8_t         nseg        = 0
    cdef int8_t         aln_flag    = 0
    cdef int32_t        nmatch      = 0
    cdef int32_t        qpos        = 0
    cdef int32_t        rpos, qend, idx, op, l, k, n, m
    cdef uint32_t       *cigar_p

    n = src.core.n_cigar
    m = n-1
    rpos    = src.core.pos
    cigar_p = bam_get_cigar(src)
    # traverse alignment's CIGAR
    for k in range(n):
        op = cigar_p[k] & BAM_CIGAR_MASK
        l  = cigar_p[k] >> BAM_CIGAR_SHIFT
        if is_match(op):
            qpos   += l
            rpos   += l
            nmatch += l
        elif is_clip_or_insert(op):
            if l >= minl:
                idx  = N + nseg
                qend = qpos + l
                if k==0:
                    segs[idx].qst      = qpos
                    segs[idx].qed      = qend
                    segs[idx].rpos     = rpos
                    segs[idx].sflag    = LEFT_CLIP
                    segs[idx].ith      = ith
                    segs[idx].overhang = nmatch
                    aln_flag |= LEFT_CLIP
                elif k==m:
                    segs[idx].qst      = qpos
                    segs[idx].qed      = qend
                    segs[idx].rpos     = rpos
                    segs[idx].sflag    = RIGHT_CLIP
                    segs[idx].ith      = ith
                    segs[idx].overhang = nmatch
                    aln_flag |= RIGHT_CLIP
                else:
                    segs[idx].qst      = qpos
                    segs[idx].qed      = qend
                    segs[idx].rpos     = rpos
                    segs[idx].sflag    = MID_INSERT
                    segs[idx].ith      = ith
                    segs[idx].overhang = nmatch
                    aln_flag |= MID_INSERT
                ith  += 1
                nseg += 1
            qpos += l
        elif is_del_or_skip(op):
            rpos += l
    
    if nseg > 0:
        idx  = N+nseg
        for k from N <= k < idx:
            segs[k].flag    = src.core.flag
            segs[k].mapq    = src.core.qual
            segs[k].rflag   = aln_flag
            segs[k].offset  = offset
            segs[k].refst   = src.core.pos
            segs[k].refed   = rpos
            segs[k].nseg    = nseg
            segs[k].nmatch  = nmatch
    
    return nseg
#
# ---------------------------------------------------------------
#
#####################
### Previous code ###
#####################

# cdef inline char* getSequenceInRange(bam1_t   *src,
#                                      uint32_t start,
#                                      uint32_t end):
#     '''return python string of the sequence in a bam1_t object.'''
#     cdef uint8_t *p
#     cdef uint32_t k
#     cdef char *s
#     cdef bytes seq
    
#     seq = PyBytes_FromStringAndSize(NULL, end - start)
#     s   = <char*>seq
#     p   = bam_get_seq(src)

#     for k from start <= k < end:
#         # equivalent to seq_nt16_str[bam1_seqi(s, i)] (see bam.c)
#         # note: do not use string literal as it will be a python string
#         s[k-start] = seq_nt16_str[p[k//2] >> 4 * (1 - k%2) & 0xf]

#     return seq
# #
# # ---------------------------------------------------------------
# #
# cdef class InsertSegment:
#     def __init__(self):
#         # see bam_init1
#         self._delegate = <bam1_t*>calloc(1, sizeof(bam1_t))
#         if self._delegate == NULL:
#             raise MemoryError("could not allocated memory of "
#                               "{} bytes".format(sizeof(bam1_t)))
#         # allocate some memory. If size is 0, calloc does not return a
#         # pointer that can be passed to free() so allocate 40 bytes
#         # for a new read
#         self._delegate.m_data = 40
#         self._delegate.data = <uint8_t *>calloc(self._delegate.m_data, 1)
#         if self._delegate.data == NULL:
#             raise MemoryError("could not allocate memory of "
#                               "{} bytes".format(self._delegate.m_data))
#         self._delegate.l_data = 0
#         # set some data to make read approximately legit.
#         # Note, SAM writing fails with q_name of length 0
#         self._delegate.core.l_qname = 0
#         self._delegate.core.tid = -1
#         self._delegate.core.pos = -1
#         self._delegate.core.mtid = -1
#         self._delegate.core.mpos = -1

#     def __dealloc__(self):
#         bam_destroy1(self._delegate)
    
#     def __lt__(self, InsertSegment other):
#         return self.rpos < other.rpos
    
#     property qname:
#         '''the query template name (None if not present)'''
#         def __get__(self):
#             cdef bam1_t* src = self._delegate
#             cdef char* s
#             # cdef str ret
#             if src.core.l_qname == 0:
#                 return None
#             s = bam_get_qname(src)
#             # return charptr_to_str(<char *>bam_get_qname(src))
#             return s.decode(TEXT_ENCODING, ERROR_HANDLER)
    
#     property flag:
#         '''properties flag'''
#         def __get__(self):
#             return self._delegate.core.flag
    
#     property tid:
#         '''tid, the index of the reference sequence'''
#         def __get__(self):
#             return self._delegate.core.tid
    
#     property ref_start:
#         '''0-based leftmost coordinate of the alignment'''
#         def __get__(self):
#             return self._delegate.core.pos
    
#     property mapq:
#         '''mapping quality'''
#         def __get__(self):
#             cdef bam1_t * src = self._delegate
#             return bam_get_qual(src)
    
#     property q_len:
#         '''
#         the length of the query/read, corresponds to the length of 
#         the sequence supplied in the BAM/SAM file. The length of a
#         query is 0 if there is no sequence in the BAM/SAM file.
#         '''
#         def __get__(self):
#             return self._delegate.core.l_qseq
    
#     property is_mapped:
#         '''true if read itself is mapped'''
#         def __get__(self):
#             return (self.flag & BAM_FUNMAP) == 0
    
#     property is_forward:
#         '''true if read is mapped to forward strand'''
#         def __get__(self):
#             return (self.flag & BAM_FREVERSE) == 0
    
#     property is_secondary:
#         '''true if not primary alignment'''
#         def __get__(self):
#             return (self.flag & BAM_FSECONDARY) != 0
    
#     property is_supplementary:
#         '''true if this is a supplementary alignment'''
#         def __get__(self):
#             return (self.flag & BAM_FSUPPLEMENTARY) != 0
    
#     property ref_len:
#         '''
#         aligned length of the read on the reference genome. equal to
#         `ref_end - _delegate.core.pos`. 
#         '''
#         def __get__(self):
#             cdef bam1_t * src
#             src = self._delegate
#             if (self.flag & BAM_FUNMAP) or src.core.n_cigar == 0:
#                 return None
#             return self.ref_end - self._delegate.core.pos
    
#     cpdef int64_t get_tag_i(self, str tag):
#         '''Get an integer aux value

#         retrieves integer aux value from the optional alignment sec-
#         tion.

#         Parameters:
#         -----------
#         tag: str
#             a two-letter tag denoting the field. can be 'NM', 'AS', 
#             'ms', 'cm', 's1', 's2', 'rl'
        
#         Returns:
#         --------
#         value: int64_t
#             an integer value correspond to the aux tag.
#         '''
#         cdef bam1_t  *src
#         cdef uint8_t *v
#         cdef bytes   btag
#         cdef int64_t value

#         src = self._delegate
#         btag = tag.encode(TEXT_ENCODING, ERROR_HANDLER)
#         v = bam_aux_get(src, btag)
#         if v == NULL:
#             raise KeyError("tag '%s' not present" % tag)
        
#         value = bam_aux2i(v)
#         return value
    
#     cpdef double get_tag_f(self, str tag):
#         '''Get an double aux value

#         retrieves double aux value from the optional alignment sec-
#         tion.

#         Parameters:
#         -----------
#         tag: str
#             a two-letter tag denoting the field. can be 'de'
        
#         Returns:
#         --------
#         value: double
#             an double value correspond to the aux tag.
#         '''
#         cdef bam1_t  *src
#         cdef uint8_t *v
#         cdef bytes   btag
#         cdef double  value

#         src = self._delegate
#         btag = tag.encode(TEXT_ENCODING, ERROR_HANDLER)
#         v = bam_aux_get(src, btag)
#         if v == NULL:
#             raise KeyError("tag '%s' not present" % tag)
        
#         value = bam_aux2f(v)
#         return value
    
#     cpdef str get_seq(self, int start=-1, int end=-1):
#         '''Get query sequence in specified region'''
#         cdef bam1_t *src = self._delegate
#         cdef int lqseq   = src.core.l_qseq
#         cdef char* s
        
#         if lqseq == 0:
#             return None
#         if start < 0:
#             start = 0
#         if end < 0:
#             end = lqseq

#         # s = charptr_to_str(getSequenceInRange(src, start, end))
#         s = getSequenceInRange(src, start, end)
#         return s.decode(TEXT_ENCODING, ERROR_HANDLER)
    
#     cpdef tuple trim(self, int tsize):
#         cdef int t_qstart, t_qend

#         t_qstart = self.qstart - tsize
#         t_qend   = self.qend + tsize

#         if t_qstart < 0:
#             t_qstart = 0
#         if t_qend > self.q_len:
#             t_qend = self.q_len
        
#         # q_st_t = qstart - t_qstart  # EXP: "query start" of the I/S fragment in the trimmed sequence
#         # q_en_t = t_qstart + (self.qend-self.qstart) # EXP: "query end" of the I/S fragment in the trimmed sequence
#         return t_qstart, t_qend

# #
# # ---------------------------------------------------------------
# #
# cdef InsertSegment makeInsertSegment(bam1_t *src,
#                                      int32_t qstart,
#                                      int32_t qend,
#                                      int32_t rpos,
#                                      int16_t stype):
#     '''return an InsertSegment object constructed from `src`
    
#     Parameters:
#     -----------
#         src: bam1_t
#             a bam1_t pointer, point to a single alignment record.

#         qstart: int32_t
#             query start of the insert/clip segment, include.

#         qend: int32_t
#             query end of the insert/clip segment, not include.

#         rpos: int32_t
#             reference position of the insert/clip segment, corres-
#             ponds to `qend`.
            
#         stype: int16_t
#             segment type, 0x1:left-clip, 0x2:mid-insert, 0x4:right-
#             clip
    
#     Returns:
#     --------
#         dest: InsertSegment
#             a InsertSegment object, with independent copy of `src`.
#     '''
#     # note that the following does not call __init__
#     cdef InsertSegment dest = InsertSegment.__new__(InsertSegment)
#     dest._delegate = bam_dup1(src)
#     dest.qstart = qstart
#     dest.qend = qend
#     dest.rpos = rpos
#     dest.stype = stype
#     return dest
# #
# # ---------------------------------------------------------------
# #
# cdef int parse_cigar1(bam1_t *src,
#                      list tmp_segl,
#                      int32_t minl=50):
#     '''parse CIGAR

#     traverse through alignment's CIGAR, extract clip/insert segment.

#     Parameters:
#     -----------
#         src: bam1_t
#             single record pointer read from BAM file.

#         tmp_segl: list
#             temporary list used to stored the new InsertSegment ins
#             tances.
        
#         minl: int32_t
#             minimum segment length, S/I CIGAR option with length >=
#             minl will be extracted.
    
#     Returns:
#     --------
#         nseg: int16_t
#             number of insert segments extracted from this alignment
#     '''
#     cdef int16_t        nseg    = 0
#     cdef int16_t        otype   = 0
#     cdef int32_t        qpos    = 0
#     cdef int32_t        rpos, qend
#     cdef int32_t        op, l, k, n, m
#     cdef uint32_t       *cigar_p
#     cdef InsertSegment  iseg

#     n = src.core.n_cigar
#     m = n-1
#     cigar_p = bam_get_cigar(src)
#     rpos = src.core.pos
#     # traverse alignment's CIGAR
#     for k from 0 <= k < n:
#         op = cigar_p[k] & BAM_CIGAR_MASK
#         l = cigar_p[k] >> BAM_CIGAR_SHIFT
#         if is_match(op):
#             qpos += l
#             rpos += l
#         elif is_clip_or_insert(op):
#             if l >= minl:
#                 qend = qpos+l
#                 if k==0:
#                     iseg   = makeInsertSegment(src, qpos, qend, rpos, LEFT_CLIP)
#                     otype |= LEFT_CLIP
#                 elif k==m:
#                     iseg   = makeInsertSegment(src, qpos, qend, rpos, RIGHT_CLIP)
#                     otype |= RIGHT_CLIP
#                 else:
#                     iseg   = makeInsertSegment(src, qpos, qend, rpos, MID_INSERT)
#                     otype |= MID_INSERT
#                 tmp_segl.append(iseg)
#                 nseg += 1
#             qpos += l
#         elif is_del_or_skip(op):
#             rpos += l
    
#     if nseg == 1:
#         iseg = tmp_segl[-1]
#         compute_feature_single(cigar_p, n, iseg, rpos, otype, nseg)
#     elif nseg > 1:
#         compute_feature_multi(cigar_p, n, tmp_segl, rpos, otype, nseg)
    
#     return nseg
# #
# # ---------------------------------------------------------------
# #
# cdef void compute_feature_single(uint32_t *cigar_p,
#                                  uint32_t n,
#                                  InsertSegment seg,
#                                  int32_t  rpos,
#                                  int16_t  otype,
#                                  int16_t  nseg):
#     '''compute features for single Insertsegment
    
#     compute features for single Insertsegment, including:
#     1. stype:    16-bit integer, keep segment-type(the lower 3 bits),
#                  alignment-type(lower 4-6 bits), number-of-segments(
#                  higher 10 bits).

#     2. ref_end:  the rightmost base position of an alignment on the 
#                  reference genome (the coordinate of the first base 
#                  after the alignment, exclusive, 0-based).

#     3. q_start:  query start of the alignment (inclusive, 0-based.)

#     4. q_end:    query end of the alignment (exclusive, 0-based.)

#     5. overhang: the anchor length of the segment. (for clip-type s-
#                  egment, equal to reference-alignment-length. for i-
#                  nsert-type segment, choose the smaller one.)

#     Parameters:
#     -----------
#         cigar_p: uint32_t *
#             pointer of the alignment's cigar data.

#         n: uint32_t
#             number of cigar.

#         seg: InsertSegment
#             InsertSegment object.

#         rpos: int32_t
#             the last reference position after traversing all the cigar.

#         otype: int16_t
#             type of the alignment, just like SAM flags. (left-clip:0x1,
#             mid-insert:0x2, right-clip:0x4).

#         nseg: int16_t
#             number of segment extracted from the same alignment
#     '''
#     cdef int32_t q_start, q_end
#     cdef int32_t op, l, overhang, overhang1
#     cdef int32_t pos      = seg._delegate.core.pos
#     cdef int16_t stype    = seg.stype
#     cdef int32_t seg_rpos = seg.rpos

#     # compute q_start
#     q_start = 0
#     op = cigar_p[0] & BAM_CIGAR_MASK
#     if op==4:
#         l = cigar_p[0] >> BAM_CIGAR_SHIFT
#         q_start = l
    
#     # compute q_end
#     q_end = seg._delegate.core.l_qseq
#     op = cigar_p[n-1] & BAM_CIGAR_MASK
#     if op==4:
#         l = cigar_p[n-1] >> BAM_CIGAR_SHIFT
#         q_end -= l
    
#     # compute overhang
#     overhang = rpos - pos
#     if stype & MID_INSERT:
#         overhang = rpos - seg_rpos
#         overhang1 = seg_rpos - pos
#         if overhang > overhang1:
#             overhang = overhang1

#     # fill features of segment
#     seg.stype    = (nseg << 6) | (otype << 3) | stype
#     seg.ref_end  = rpos
#     seg.q_start  = q_start
#     seg.q_end    = q_end
#     seg.overhang = overhang
# #
# # ---------------------------------------------------------------
# #
# cdef void compute_feature_multi(uint32_t *cigar_p,
#                                 uint32_t n,
#                                 list     tmp_segl,
#                                 int32_t  rpos,
#                                 int16_t  otype,
#                                 int16_t  nseg):
#     '''compute features for multiple Insertsegment

#     compute features for multiple Insertsegment, including:
#     1. stype:    16-bit integer, keep segment-type(the lower 3 bits),
#                  alignment-type(lower 4-6 bits), number-of-segments(
#                  higher 10 bits).

#     2. ref_end:  the rightmost base position of an alignment on the 
#                  reference genome (the coordinate of the first base 
#                  after the alignment, exclusive, 0-based).

#     3. q_start:  query start of the alignment (inclusive, 0-based.)

#     4. q_end:    query end of the alignment (exclusive, 0-based.)

#     5. overhang: the anchor length of the segment. (for clip-type s-
#                  egment, equal to reference-alignment-length. for i-
#                  nsert-type segment, choose the smaller one.)

#     6. ldist:    reference distance to the upstream segment from the
#                  same alignment.

#     7. rdist:    reference distance to the downstream segment.

#     Parameters:
#     -----------
#         cigar_p: uint32_t *
#             pointer of the alignment's cigar data.

#         n: uint32_t
#             number of cigar.

#         tmp_segl: list
#             list of InsertSegment objects.

#         rpos: int32_t
#             the last reference position after traversing all the cigar.

#         otype: int16_t
#             type of the alignment, just like SAM flags. (left-clip:0x1,
#             mid-insert:0x2, right-clip:0x4).
            
#         nseg: int16_t
#             number of segment extracted from the same alignment
#     '''
#     cdef InsertSegment  seg = tmp_segl[-nseg], p_seg
#     cdef int32_t    pos     = seg._delegate.core.pos
#     cdef int16_t    stype   = seg.stype
#     cdef int32_t    s_rpos  = seg.rpos, p_rpos
#     cdef int32_t    op, l, overhang, overhang1
#     cdef int32_t    q_start, q_end, dist

#     # compute q_start
#     q_start = 0
#     op = cigar_p[0] & BAM_CIGAR_MASK
#     if op==4:
#         l = cigar_p[0] >> BAM_CIGAR_SHIFT
#         q_start = l
    
#     # compute q_end
#     q_end = seg._delegate.core.l_qseq
#     op = cigar_p[n-1] & BAM_CIGAR_MASK
#     if op==4:
#         l = cigar_p[n-1] >> BAM_CIGAR_SHIFT
#         q_end = q_end - l
    
#     # compute overhang for the first segment
#     overhang =  rpos - pos
#     if stype & MID_INSERT:
#         overhang  = rpos - s_rpos
#         overhang1 = s_rpos - pos
#         if overhang > overhang1:
#             overhang = overhang1
#     # fill features of the first segment
#     seg.stype    = (nseg << 6) | (otype << 3) | stype
#     seg.ref_end  = rpos
#     seg.q_start  = q_start
#     seg.q_end    = q_end
#     seg.overhang = overhang


#     # compute overhang & distance
#     cdef int32_t i
#     for i in range(-nseg+1, 0):
#         # compute overhang
#         seg      =  tmp_segl[i]
#         s_rpos   =  seg.rpos
#         stype    =  seg.stype
#         overhang =  rpos - pos
#         if stype & MID_INSERT:
#             overhang  = rpos - s_rpos
#             overhang1 = s_rpos - pos
#             if overhang > overhang1:
#                 overhang = overhang1

#         # compute distance
#         p_seg   = tmp_segl[i-1]
#         p_rpos  = p_seg.rpos
#         dist    = s_rpos - p_rpos
        
#         # fill features of segment
#         seg.stype    = (nseg << 6) | (otype << 3) | stype
#         seg.ref_end  = rpos
#         seg.q_start  = q_start
#         seg.q_end    = q_end
#         seg.overhang = overhang
#         seg.ldist    = dist
#         p_seg.rdist  = dist