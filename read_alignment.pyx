from pysam import AlignmentFile
from uuid import uuid4
from collections import OrderedDict, defaultdict
from bx.intervals.intersection import Intersecter, Interval
import mappy as mp
from subprocess import Popen, PIPE, STDOUT, DEVNULL
import os



##
cdef class Segment:
    cdef public:
        str segname, te_name
        int q_st, q_en
        int q_st_t, q_en_t
        int qlen, seg_len
        object r_st, r_en
        int rpos, mapq
        int te_mapq, te_r_st, te_r_en, te_q_st, te_q_en, te_orient
        object te_aln, read
        bint is_reverse

    def __init__(self, object read, str segname, int q_st, int q_en, object r_st, object r_en, int rpos, int q_st_t, int q_en_t, int qlen):
        self.segname = segname
        self.q_st = q_st # 0-based, include. Left-most position of the insert/clip segment on the query sequence.
        self.q_en = q_en # 0-based, exclude. Left-most position of the insert/clip segment on the query sequence.
        self.qlen = qlen
        self.seg_len = q_en - q_st

        self.r_st = r_st # reference segment start. reference position of insert/clip segment start.
        self.r_en = r_en # reference segment end. reference position of insert/clip segment end.
        self.rpos = rpos
        self.mapq = read.mapping_quality
        self.is_reverse = read.is_reverse

        self.read = read


##
cdef class Cluster:
    cdef public:
        str c_id, consensus, tsd, ref_name
        list tes, segs, bp

    def __init__(self, str c_id, str ref_name):
        self.c_id = c_id
        self.ref_name = ref_name
        self.segs = []
        self.tes = []

    cpdef add_seg(self, Segment seg):
        """ 向Cluster对象的segs当中添加周围的Segment """
        self.segs.append(seg)

    cpdef find_bp(self):
        """ find break points """
        cdef:
            list bps
            Segment s
        
        bps = sorted([s.rpos for s in self.segs])
        self.bp = [bps[0], bps[-1]]
    
    cpdef str cluster2bed(self):
        return '{}\t{}\t{}\n'.format(self.ref_name, self.bp[0], self.bp[-1])
    
    # cpdef map_te(self, str te_idx, str p="splice", int n=5):
    #     """ align trimmed segment sequence(s) to transposon consensus sequence"""
    #     cdef:
    #         object aln=mp.Aligner(te_idx, preset=p, best_n=n), hit
    #         Segment seg
        
    #     for seg in self.segs:
    #         if seg.qseq_t:
    #             for hit in aln.map(seg.qseq_t):
    #                 self.tes.append((seg.qname, hit))

    def make_consensus(self):
        """ define consensus sequence """
        pass

    def find_tsd(self):
        """ find TSD sequence """
        pass


##
cdef trim(int qlen, str seq, int q_st, int q_en, int op_len, int flanksize):
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

    return seq_t, q_st_t, q_en_t


##
cdef extract_seg(object read, object seg_fa, int flanksize=200, unsigned int min_len=200, unsigned int max_len=10000):
    """ 根据输入的alignment (read)，解析其CIGAR字段并提取alignment当中的clip或
        insert segment。最后返回一个list (segs)，其中的每个元素都是Segment对象 """
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
                    segname = '{}_{}_{}_{}_{}_r'.format(read.query_name, read.reference_name, rpos, q_st, q_en)
                    # segname = read.qname + '_' + read.reference_name + '_r'
                else:
                    segname = '{}_{}_{}_{}_{}_l'.format(read.query_name, read.reference_name, rpos, q_st, q_en)
            elif op == 1: # I
                segname = '{}_{}_{}_{}_{}_m'.format(read.query_name, read.reference_name, rpos, q_st, q_en)

            # trim the query sequence and define the trimmed query start & end of the segment
            if seq:
                seq_t, q_st_t, q_en_t = trim(qlen, seq, q_st, q_en, op_len, flanksize)
                
                # write out trimmed sequence
                seg = Segment(read, segname, q_st, q_en, r_st, r_en, rpos, q_st_t, q_en_t, qlen)
                seg_fa.write('>{}\n{}\n'.format(segname, seq_t))
            else:
                seg = Segment(read, segname, q_st, q_en, r_st, r_en, rpos, 0, 0, qlen)

            segs.append(seg)
            
        if op in (0, 1, 4, 7, 8):
            start_idx += op_len
    
    return segs # list of Segment


##
cpdef build_cluster(str bam, str chrom, int flank=200):
    cdef:
        object bam_file = AlignmentFile(bam, 'rb'), seg_fa = open(chrom+".tmp.fa", "w")
        dict id2cluster = {}, id2cluster_p
        list segs, c_ids
        Segment seg
        Cluster cluster
        str c_id
        object tree = Intersecter(), read
    
    for read in bam_file.fetch(chrom):
        segs = extract_seg(read, seg_fa)
        # if read.query_sequence is None or read.qual is None:
        #     print(read.query_name, read.is_secondary, read.is_supplementary)

        for seg in segs:
            c_ids = [x.value for x in tree.find(seg.rpos-1, seg.rpos+1)]
            if c_ids:
                for c_id in c_ids:
                    id2cluster[c_id].add_seg(seg)
            else:        
                c_id = str(uuid4())
                cluster = Cluster(c_id, chrom)
                id2cluster[c_id] = cluster
                id2cluster[c_id].add_seg(seg)
                tree.add_interval(Interval(seg.rpos-flank, seg.rpos+flank, value=c_id))
    
    seg_fa.close()
    return id2cluster


##
cpdef process_cluster(dict id2cluster, str chrom, str te_idx):
    cdef:
        str c_id, te_aln_path
        Cluster cluster
        object fout = open("{}.tmp.bed".format(chrom), "w")
    
    te_aln_path = align_mm2(te_idx, chrom + ".tmp.fa", thread_mm2=10)
    for c_id in id2cluster:
        id2cluster[c_id].find_bp()
        # id2cluster[c_id].map_te(te_idx)
        fout.write(id2cluster[c_id].cluster2bed())
    
    fout.close()
    return id2cluster

##
cdef align_mm2(str ref, str query, int thread_mm2, str preset="map-ont"):
    cdef:
        int exitcode
        str wk_dir = os.path.dirname(query)
        str q_p = os.path.basename(query).rsplit('.', 1)[0]
        str ref_p = os.path.basename(ref).rsplit('.', 1)[0]
        str fout = os.path.join(wk_dir, ref_p + '_' + q_p + '.bam')
        list cmd_mm2 = ['minimap2' + ' -t ' + str(thread_mm2) + ' -ax ' + preset + ' ' + ref + ' ' + query + ' | ' + 'samtools view -Shb - -o ' + fout]
    
    align_proc = Popen(cmd_mm2, stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = align_proc.wait()
    if exitcode != 0:
        raise Exception("Error: minimap2 alignment for {} to {} failed".format(q_p, ref_p))
    
    return fout
