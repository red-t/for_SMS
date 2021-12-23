from pysam import AlignmentFile
from uuid import uuid4
from collections import OrderedDict, defaultdict
from bx.intervals.intersection import Intersecter, Interval


##
cdef class Segment:
    cdef public:
        str qname, qseq, qual, ref_name, qseq_t, qual_t, te_name, segname
        int qsstart, qsend, qlen, seg_len, rspos, mapq, qsstart_t, qsend_t, te_mapq, te_rstart, te_rend, te_qstart, te_orient
        object te_aln, read, rsstart, rsend
        bint is_reverse

    def __init__(self, object read, int q_ss, int q_se, object r_ss, object r_se, int rspos):
        self.qname = read.query_name
        self.qseq = read.query_sequence
        self.qual = read.qual
        self.qsstart = q_ss # 0-based, include. Left-most position of the insert/clip segment on the query sequence.
        self.qsend = q_se # 0-based, exclude. Left-most position of the insert/clip segment on the query sequence.
        self.qlen = read.infer_query_length()
        self.seg_len = self.qsend - self.qsstart

        self.ref_name = read.reference_name
        self.rsstart = r_ss # reference segment start. reference position of insert/clip segment start.
        self.rsend = r_se # reference segment end. reference position of insert/clip segment end.
        self.rspos = rspos
        self.mapq = read.mapping_quality
        self.is_reverse = read.is_reverse

        self.read = read
    
    cdef trim(self, int flanksize=200):
        """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。 """
        cdef:
            int t_start, t_end

        t_start = self.qsstart - flanksize
        t_end = self.qsend + flanksize
        
        if t_start < 0:
            t_start = 0
        if t_end > self.qlen:
            t_end = self.qlen
        
        self.qsstart_t = self.qsstart - t_start
        self.qsend_t = self.qsstart_t + self.seg_len
        if self.qseq:
            self.qseq_t = self.qseq[t_start:t_end]
        
        if self.qual:
            self.qual_t = self.qual[t_start:t_end]
    
    cdef str seg2fq(self):
        return '@{}\n{}\n+\n{}\n'.format(self.qname, self.qseq_t, self.qual_t)
    
    cdef str seg2fa(self):
        return '>{}_{}\n{}\n'.format(self.qname, self.op, self.qseq_t)
    
    cdef parse_type(self, int op):
        if op == 4: # S
            if self.qsstart == 0:
                self.segname = self.qname + '~' + str(uuid4()) + '~r'
            else:
                self.segname = self.qname + '~' + str(uuid4()) + '~l'
        elif op == 1: # I
            self.segname = self.qname + str(uuid4()) + '~m'


##
cdef class Cluster:
    cdef public:
        str c_id, consensus, tsd, ref_name
        list tes, segs, bp

    def __init__(self, str c_id, str ref_name):
        self.c_id = c_id
        self.ref_name = ref_name
        self.segs = []

    cpdef add_seg(self, Segment seg):
        """ 向Cluster对象的segs当中添加周围的Segment """
        self.segs.append(seg)

    cpdef find_bp(self):
        """ find break points """
        cdef:
            list bps
            Segment s
        
        bps = sorted([s.rspos for s in self.segs])
        self.bp = [bps[0], bps[-1]]
    
    cpdef str cluster2bed(self):
        return '{}\t{}\t{}\n'.format(self.ref_name, self.bp[0], self.bp[-1])

    def make_consensus(self):
        """ define consensus sequence """
        pass

    def find_tsd(self):
        """ find TSD sequence """
        pass


##
cdef extract_seg(object read, unsigned int min_len=200, unsigned int max_len=10000):
    """ 根据输入的alignment (read)，解析其CIGAR字段并提取alignment当中的clip或
        insert segment。最后返回一个list (segs)，其中的每个元素都是Segment对象 """
    cdef:
        int start_idx=0, op, length, q_ss, q_se
        object r_ss, r_se
        list segs=[]
        dict ap=dict(read.get_aligned_pairs())
        Segment seg

    ap[-1] = None
    ap[read.infer_query_length()] = None

    for op, length in read.cigartuples: # op:CIGAR operation; length: length of this operation
        if op in (1, 4) and length >= min_len: # I or S
            if length > max_len:
                continue

            q_ss = start_idx # 0-based, include, query insert start index of each I/S fragment.
            q_se = start_idx + length # 0-based, exclude, query insert end index of each I/S fragment.

            r_ss = ap[q_ss - 1]
            r_se = ap[q_se]

            if r_ss:
                rspos = r_ss
            elif r_se:
                rspos = r_se
            else:
                continue
                
            seg = Segment(read, q_ss, q_se, r_ss, r_se, rspos)
            seg.trim()
            seg.parse_type(op)

            segs.append(seg)
            
        if op in (0, 1, 4, 7, 8):
            start_idx += length
    
    return segs # list of Segment


##
cpdef build_cluster(str bam, str chrom, int flank=200):
    cdef:
        object bam_file = AlignmentFile(bam, 'rb')
        dict id2cluster = {}, id2cluster_p
        list segs, c_ids
        Segment seg
        Cluster cluster
        str c_id
        object tree = Intersecter(), read
    
    for read in bam_file.fetch(chrom):
        segs = extract_seg(read)
        # if read.query_sequence is None or read.qual is None:
        #     print(read.query_name, read.is_secondary, read.is_supplementary)

        for seg in segs:
            c_ids = [x.value for x in tree.find(seg.rspos-1, seg.rspos+1)]
            if c_ids:
                for c_id in c_ids:
                    id2cluster[c_id].add_seg(seg)
            else:        
                c_id = str(uuid4())
                cluster = Cluster(c_id, seg.ref_name)
                id2cluster[c_id] = cluster
                id2cluster[c_id].add_seg(seg)
                tree.add_interval(Interval(seg.rspos-flank, seg.rspos+flank, value=c_id))

    return id2cluster


##
cpdef process_cluster(dict id2cluster, str chrom):
    cdef:
        str c_id
        Cluster cluster
        object fout = open("{}.tmp.bed".format(chrom), "w")
    
    for c_id in id2cluster:
        id2cluster[c_id].find_bp()
        fout.write(id2cluster[c_id].cluster2bed())
    
    fout.close()
    return id2cluster


# cdef trim(int qlen, str seq, int q_st, int q_en, int op_len, int flanksize):
#     """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。 """
#     cdef:
#         int st_t, en_t
#         int q_st_t, q_en_t
#         str seq_t

#     st_t = q_st - flanksize
#     en_t = q_en + flanksize

#     if st_t < 0:
#         st_t = 0
#     if en_t > qlen:
#         en_t = qlen

#     q_st_t = q_st - st_t
#     q_en_t = q_st_t + op_len
#     seq_t = seq[st_t:en_t]

#     return seq_t, q_st_t, q_en_t