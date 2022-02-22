from pysam import AlignmentFile
from uuid import uuid4
from collections import OrderedDict, defaultdict, Counter
from bx.intervals.intersection import Intersecter, Interval
from mappy import revcomp
from subprocess import Popen, PIPE, STDOUT, DEVNULL
import os


## Segment ##
cdef class Segment:
    cdef public:
        str segname, seq_t
        int strand
        int q_st, q_en
        int q_st_t, q_en_t
        object r_st, r_en
        int rpos
        list type_list, cord_list, orients
        int n_aln
        object read

    def __init__(self, object read, str segname, bint strand, int q_st, int q_en, object r_st, object r_en, int rpos, int q_st_t, int q_en_t, str seq_t=''):
        self.segname = segname
        self.strand = int(strand)

        self.q_st = q_st # 0-based, include. Left-most position of the insert/clip segment on the query sequence.
        self.q_en = q_en # 0-based, exclude. Right-most position of the insert/clip segment on the query sequence.

        self.seq_t = seq_t
        self.q_st_t = q_st_t
        self.q_en_t = q_en_t

        self.r_st = r_st # reference segment start. reference position of insert/clip segment start.
        self.r_en = r_en # reference segment end. reference position of insert/clip segment end.
        self.rpos = rpos
        # self.mapq = read.mapping_quality
        self.type_list = []
        self.cord_list = []
        self.orients = []
        self.read = read
    
    cpdef add_te(self, str te_n, int q_st, int q_en, int r_st, int r_en, int strand):
        self.type_list.append(te_n)
        self.cord_list.append((q_st, q_en, r_st, r_en, strand))
        self.orients.append(abs(self.strand - strand))
        self.n_aln += 1
    
    @property
    def n_type(self):
        return len(set(self.type_list))



## Cluster ##
cdef class Cluster:
    cdef public:
        int orient
        str c_id, consensus, tsd, ref_name, te_type
        list type_list, seg_list, orients, bp

    def __init__(self, str c_id, str ref_name):
        self.c_id = c_id
        self.ref_name = ref_name
        self.seg_list = []
        self.type_list = []
        self.orients = []

    cpdef add_seg(self, Segment seg):
        """ 向Cluster对象的segs当中添加周围的Segment """
        self.seg_list.append(seg)
        self.type_list.extend(seg.type_list)
        self.orients.extend(seg.orients)

    cpdef find_bp(self):
        """ find break points """
        cdef:
            list bps
            Segment s
        
        bps = sorted([s.rpos for s in self.seg_list])
        self.bp = [bps[0], bps[-1]]

    cpdef te_stat(self):
        self.orient = Counter(self.orients).most_common()[0][0]
        self.te_type = Counter(self.type_list).most_common()[0][0]


    cpdef str cluster2bed(self):
        cdef:
            str strand

        if self.orient:
            strand = '-'
        else:
            strand = '+'

        return '{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.ref_name, self.bp[0]+1, self.bp[-1]+1, self.te_type, 0, strand)

    def make_consensus(self):
        """ define consensus sequence """
        pass

    def find_tsd(self):
        """ find TSD sequence """
        pass


## trim ##
cdef tuple trim(bint is_reverse, int qlen, str seq, int q_st, int q_en, int op_len, int flanksize):
    """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。 """
    cdef:
        int st_t, en_t
        int q_st_t, q_en_t
        str seq_t

    st_t = q_st - flanksize
    en_t = q_en + flanksize

    if st_t < 0:
        st_t = 0
    if en_t > qlen:
        en_t = qlen

    q_st_t = q_st - st_t
    q_en_t = q_st_t + op_len
    seq_t = seq[st_t:en_t]
    if is_reverse:
        seq_t = revcomp(seq_t)

    return seq_t, q_st_t, q_en_t


## extract_seg ##
cdef list extract_seg(object read, object seg_fa, int flanksize=200, unsigned int min_len=200, unsigned int max_len=10000):
    """ 根据输入的alignment，解析其CIGAR字段并提取alignment当中的clip或
        insert segment。最后返回list of tuple, (segname, seg) """
    cdef:
        int start_idx=0
        int op, op_len # `op`:CIGAR operation; `op_len`: length of this operation
        int q_st, q_en # 0-based, query start/end index of each I/S segment. `q_st` include, `q_en` exclude
        int q_st_t, q_en_t # trimmed `q_st` & `q_en`
        object r_st, r_en # 0-based, reference start/end index of each I/S segment. aligned pairs of `q_st-1` and `q_en` respectively
        int qlen = read.infer_query_length() # len of the query sequence (read)
        list segs=[]
        dict ap=dict(read.get_aligned_pairs())
        Segment seg
        str segname # name of the I/S segment
        str seq_t, seq = read.query_sequence

    ap[-1] = None
    ap[qlen] = None

    for op, op_len in read.cigartuples:
        if op in (1, 4) and op_len >= min_len: # I or S
            if op_len > max_len:
                continue
            
            # find query start & end of the segment
            q_st = start_idx 
            q_en = start_idx + op_len

            # find reference start & end of the segment
            r_st = ap[q_st - 1]
            r_en = ap[q_en]

            # define reference position of the segment
            if r_st:
                rpos = r_st
            elif r_en:
                rpos = r_en
            else:
                continue
            
            # generate name of the segment
            if op == 4: # S
                if q_st == 0:
                    segname = '{}_{}_{}_{}_{}_{}_l'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)
                else:
                    segname = '{}_{}_{}_{}_{}_{}_r'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)
            elif op == 1: # I
                segname = '{}_{}_{}_{}_{}_{}_m'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)

            # trim the query sequence and define the trimmed query start & end of the segment
            if seq:
                seq_t, q_st_t, q_en_t = trim(read.is_reverse, qlen, seq, q_st, q_en, op_len, flanksize)
                
                # write out trimmed sequence
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, q_st_t, q_en_t, seq_t)
                seg_fa.write('>{}\n{}\n'.format(segname, seq_t))
            else:
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, 0, 0)

            segs.append((segname, seg))
            
        if op in (0, 1, 4, 7, 8):
            start_idx += op_len

    return segs # list of tuples, (segname, seg)


## align_mm2 ##
cdef str align_mm2(str ref, str query, int thread_mm2, str preset="map-ont"):
    cdef:
        int exitcode
        str wk_dir = os.path.dirname(query)
        str q_p = os.path.basename(query).rsplit('.', 1)[0]
        str ref_p = os.path.basename(ref).rsplit('.', 1)[0]
        str out_path = os.path.join(wk_dir, ref_p + '_' + q_p + '.bam')
        list cmd_mm2 = ['minimap2', '-t', str(thread_mm2), '-aYx', preset, '--secondary=no', ref, query, '|', 'samtools', 'view', '-bhSG 4', '-', '|', 'samtools', 'sort', '-o', out_path, '-']
        list cmd_idx = ['samtools index ' + out_path]
    
    align_proc = Popen([" ".join(cmd_mm2)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = align_proc.wait()
    if exitcode != 0:
        raise Exception("Error: minimap2 alignment for {} to {} failed".format(q_p, ref_p))
    
    idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = idx_proc.wait()
    if exitcode != 0:
        raise Exception("Error: samtools index for {} failed".format(out_path))
    return out_path

## parse_te ##
cdef parse_te_aln(dict seg_dict, str te_bam):
    cdef:
        object te_aln = AlignmentFile(te_bam, "rb"), r
        str seg_n
        int q_st, q_en
    
    for r in te_aln:
        if not r.is_unmapped:
            
            seg_n = r.query_name
            if r.is_reverse:
                q_st = r.query_length - r.qend
                q_en = q_st + r.qlen           
            else:
                q_st = r.qstart
                q_en = r.qend

            seg_dict[seg_n].add_te(r.reference_name, q_st, q_en, r.reference_start, r.reference_end, int(r.is_reverse))

## collect_seg ##
cdef dict collect_seg(str ref_bam, str chrom, str te_idx, int flank=200):
    cdef:
        object ref_aln = AlignmentFile(ref_bam, 'rb'), seg_fa = open(chrom+".tmp.fa", "w")
        dict seg_dict = {}
        list segs = []
        object read
        str te_bam
    
    # extract segment from reference genome alignment
    for read in ref_aln.fetch(chrom):
        if not read.is_secondary:
            segs.extend(extract_seg(read, seg_fa))

    seg_fa.close()
    seg_dict = dict(segs) # should we remove segs after this?

    # align segment sequence to transposon consensus sequence
    te_bam = align_mm2(te_idx, chrom + ".tmp.fa", thread_mm2=10)

    # parse TE alignment of segment sequence
    parse_te_aln(seg_dict, te_bam)

    return seg_dict


## build_cluster ##
cpdef dict build_cluster(str ref_bam, str chrom, str te_idx, int flank=200):
    cdef:
        str c_id
        list c_ids
        dict seg_dict
        dict cluster_dict = {}
        object tree = Intersecter()
        Segment seg
        Cluster cluster
    
    seg_dict = collect_seg(ref_bam, chrom, te_idx, flank)

    for seg in seg_dict.values():
        if seg.n_aln:
            c_ids = [x.value for x in tree.find(seg.rpos-1, seg.rpos+1)]
            if c_ids:
                for c_id in c_ids:
                    cluster_dict[c_id].add_seg(seg)
            else:        
                c_id = str(uuid4())
                cluster = Cluster(c_id, chrom)
                cluster_dict[c_id] = cluster
                cluster_dict[c_id].add_seg(seg)
                tree.add_interval(Interval(seg.rpos-flank, seg.rpos+flank, value=c_id))
        
    return cluster_dict


## process_cluster ##
cpdef process_cluster(dict cluster_dict, str chrom):
    cdef:
        str c_id
        Cluster cluster
        object fout = open("{}.tmp.bed".format(chrom), "w")
    
    for c_id in cluster_dict:
        cluster_dict[c_id].find_bp()
        cluster_dict[c_id].te_stat()
        fout.write(cluster_dict[c_id].cluster2bed())
    
    fout.close()
    return cluster_dict
