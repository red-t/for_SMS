from pysam import AlignmentFile, AlignedSegment, FastaFile
from uuid import uuid4
from collections import OrderedDict, defaultdict, Counter
from bx.intervals.intersection import Intersecter, Interval
from mappy import revcomp
from subprocess import Popen, PIPE, STDOUT, DEVNULL
import os
from consensus import get_consensus_seq
from operator import itemgetter
import re

# 目前的筛选条件：
# 1. primary alignment & supplementary alignment
# 2. insertion 长度在200-10000之间
# 3. flanksize200

# 添加
# 4. hard clip也可能是
# 5. 有些比对到genome上的质量比较差，添加mapq
# 6. te list里并不全是一类转座子
# 7. 把insert序列信息记录下来

# 8. 给TE 比对的情况打分
# 9. 筛选supporting reads和unspporting reads 计算frequency

# 7. 区分出nested insertion
# 8. 把insert序列做consensus
# 9. consensus两端跟genome比取TSD

# 10.方向问题，最终输出以genone reference 正链为参考，
#            bam文件中sequence序列是与正链一致的，所以insert seq比对到TE的方向就是insertion的方向


## Segment ##
cdef class Segment:
    cdef public:
        str segname, seq_t
        int genome_strand
        int q_st, q_en
        int q_st_t, q_en_t
        object r_st, r_en
        int rpos
        list type_list, cord_list, orients
        int n_aln
        object read

    def __init__(self, object read, str segname, bint genome_strand, int q_st, int q_en, object r_st, object r_en, int rpos, int q_st_t, int q_en_t, str seq_t=''):
        self.segname = segname
        self.genome_strand = int(genome_strand)

        self.q_st = q_st # 0-based, include. Left-most position of the insert/clip segment on the query sequence.
        self.q_en = q_en # 0-based, exclude. Right-most position of the insert/clip segment on the query sequence.

        self.seq_t = seq_t
        self.q_st_t = q_st_t
        self.q_en_t = q_en_t

        self.r_st = r_st # reference segment start. reference position of insert/clip segment start.
        self.r_en = r_en # reference segment end. reference position of insert/clip segment end.
        self.rpos = rpos
        self.type_list = []
        self.cord_list = []
        self.orients = []
        self.read = read
    
    cpdef add_te(self, str te_n, int q_st, int q_en, int r_st, int r_en, int te_strand):
        
        self.cord_list.append( (q_st, q_en, r_st, r_en, te_strand) )
        self.orients.append( te_strand )

        # segment里每条supporting reads对应的TE 类型和 方向要保持对应，方便判断nested insertion
        self.type_list.append( te_n + "_|_" + str(te_strand))
        self.n_aln += 1
    
    @property
    def n_type(self):
        return len(set(self.type_list))



## Cluster ##
cdef class Cluster:
    cdef public:
        int orient, state, supp_reads_num
        float frequency
        str c_id, consensus, consensus_seq, tsd, ref_name, te_type, insert_seq_info
        str out_path, supp_reads, unsupp_reads, span_reads
        list type_list, seg_list, orients, bp, supp_reads_list, unsupp_reads_list, span_reads_list, segname_list
        object ref_aln
        str te_idx 


    def __init__(self, str c_id, str ref_name, str out_path, object ref_aln, str te_idx):
        self.c_id = c_id
        self.ref_name = ref_name
        self.out_path = out_path
        self.seg_list = []
        self.type_list = []
        self.orients = []
        self.supp_reads_list = []
        self.unsupp_reads_list = []
        self.span_reads_list = []
        self.consensus = ''
        self.tsd = ''
        self.insert_seq_info = ''
        self.state = 1 # insertion 的是否保留的标签，有两个值，0和1，0表示这个insertion是假的，会丢掉，1就会被留下
        self.ref_aln = ref_aln
        self.frequency = 0
        self.supp_reads_num = 0
        self.segname_list = []
        self.te_idx = te_idx

    cpdef add_seg(self, Segment seg):
        """ 向Cluster对象的segs当中添加周围的Segment """
        self.seg_list.append(seg)
        self.type_list.extend(seg.type_list)
        self.orients.extend(seg.orients)
        self.supp_reads_list.append(seg.read)
        self.segname_list.append(seg.segname)

    cpdef find_bp(self):
        """ find break points """
        cdef:
            list bps
            Segment s
        
        bps = sorted([s.rpos for s in self.seg_list])
        self.bp = [bps[0], bps[-1]]

        

    cpdef te_stat(self, repeatmasker_file):
        
        cdef:
            str te_temp
            list te_temp_list=[]
        # 暂时输出一个insertion位置上所有可能的TE type
        # 但nested insertion的方向还没处理

        # 写入文件跟rmk对比去掉repeatmasker附近的insertion

        g = open( self.out_path + self.ref_name + "/" + self.c_id + "te_type.bed", 'w')

        for te_temp in list(set(self.type_list)):
            te_temp_list.append(te_temp.split('_|_')[0])
            self.orient = int(te_temp.split('_|_')[1])

            if self.orient == 1:
                te_temp_strand = '-'
            else:
                te_temp_strand = '+'
            g.write("\t".join([self.ref_name, str(self.bp[0]), str(self.bp[1]), te_temp.split('_|_')[0], str(0), te_temp_strand ]) + '\n')

        
        self.te_type = "|".join(te_temp_list)

        
        g.close()

        rm_rmk_cmd = ['bedtools', 'intersect', '-a', repeatmasker_file, '-b', self.out_path + self.ref_name + "/" + self.c_id + "te_type.bed" ,'-wa','-wb', '|','awk', '\'{if($4==$10 && $6==$12){print $4}}\'' ]
        rmk_te_list = os.popen(" ".join(rm_rmk_cmd)).readlines()
        if len(rmk_te_list) > 0:
            rmk_te = rmk_te_list[0].strip()
            if rmk_te in te_temp_list:
                te_temp_list.remove(rmk_te)
        
        self.type_list = te_temp_list
        self.te_type = "|".join(self.type_list)
        if len(self.type_list) == 0:
            self.state = 0
        
        # print(self.te_type)


        # print(">"+str(rmk_te.readlines()[0].strip()))
        # rm_rmk_proc = Popen([" ".join(rm_rmk_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        
        # self.te_type = te_temp.split('_|_')[0]
        # self.orient = int(te_temp.split('_|_')[1])

    # MrBleem
    cpdef static(self):
        """ 筛选supporting reads和unsuporting reads，计算frequency """
        # supporting reads
        cdef:
            dict black_region = {'start':1000000000,'end':0}
            int black_read_num = 0
            list support_span_reads_list = []
            # self.state insertion 的是否保留的标签，有两个值，0和1，0表示这个insertion是假的，会丢掉，1就会被留下
        for supp_read in self.supp_reads_list:
            # 处理reads比对到genome上是两端截断的情况，这种区域定义为black region
            # 如果两端截断都超过100，而且数量比较多，那么可能这个区域是genome上的reference TE
            # 作出这种处理是因为在IGV上看出这种两端截断的结构很明显
            # chr2L:2,932,687-2,936,323
            if supp_read.cigartuples[0][0] in (4,5) and supp_read.cigartuples[-1][0] in (4, 5):
                # self.state = 0
                # print(supp_read.cigartuples[0][1], supp_read.cigartuples[-1][1])
                if supp_read.cigartuples[0][1] > 100 and supp_read.cigartuples[-1][1] > 100:
                    black_read_num = black_read_num + 1
                    if black_region['start'] > supp_read.reference_start:
                        black_region['start'] = supp_read.reference_start
                    if supp_read.reference_end > black_region['end']:
                        black_region['end'] = supp_read.reference_end
                    if black_read_num >=2:
                        self.state = 0
            else:
                support_span_reads_list.append(supp_read)
        
        # 如果supporting reads里有一条reads是跨过这条black region的，那么这个位置的insertion就会被保留下来
        for supp_span_read in support_span_reads_list:
            if black_region['start'] - supp_span_read.reference_start > 100 or supp_span_read.reference_end - black_region['end'] > 100:
                self.state = 1
                continue

        
        
        # 这里是处理suporting reads里左右两边clip reads数量比例的问题
        # 如果有一条reads是跨过insertion的，多少clip reads都无所谓
        # 如果没有跨过insertion的reads，且只有一边有比较多的（>=2）clip reads,这个insertion也会丢掉
        # 如果两边都有，则计算比例，0.1883698 是算了一下，所有只有clip reads的insertion，左右两边clip reads的比例，拟合出来一个正太分布
        # 取的是正太分布的（平均值 - 3 * 标准差  ～ 平均值 + 3 * 标准差 ）
        supp_reads_list_id = list(set(['_'.join(r.split('_')[3:5]) for r in self.segname_list ]))
        # supp_reads_list_id = list(['_'.join(r.split('_')[3:5]) for r in self.segname_list ])
        supp_reads_type = [ r[-1] for r in self.segname_list ]
        self.supp_reads = str(len(supp_reads_list_id)) +"_|_"+ "_|_".join(supp_reads_list_id) + "_|_".join(self.segname_list)
        n_m = supp_reads_type.count('m')
        n_l = supp_reads_type.count('l')
        n_r = supp_reads_type.count('r')
        # print(n_m, n_l, n_r)

        if n_m == 0 :
            if n_l == 0 :
                if n_r >= 2:
                    self.state = 0
            elif n_r == 0 :
                if n_l >= 2:
                    self.state = 0
            else:
                if  abs( 0.5 - ( float(n_l) / sum([n_l,n_r])) ) >  0.5 - 0.1883698: 
                    self.state = 0
            

        
        # unsupporting reads
        if self.bp[0] - 10 <= 0 :
            start_bp = 1
        else:
            start_bp = self.bp[0] - 10

        # fetch得到insertion position区域的所有reads
        # 如果reads跨过insertion position两端各30bp以上，而且不是supporting reads，这条reads认为是unsuporting reads
        for span_read in self.ref_aln.fetch(contig=self.ref_name, start=start_bp, stop=self.bp[-1]+10):
            #print(span_read.reference_start, self.bp[0], span_read.reference_end, self.bp[-1])
            if self.bp[0] - span_read.reference_start >= 30 and span_read.reference_end - self.bp[-1] >= 30:
                span_read_id = '_'.join(span_read.query_name.split('_')[3:5])
                if span_read_id not in supp_reads_list_id:
                    self.unsupp_reads_list.append(span_read_id)
        
        # frequency
        self.frequency =float(len(supp_reads_list_id)) / ( len(supp_reads_list_id) + len(self.unsupp_reads_list))




    cpdef str cluster2bed(self):
        cdef:
            str strand

        if self.orient:
            strand = '-'
        else:
            strand = '+'
        
        nested = ''
        if len(list(set(self.type_list))) > 1:
            nested = 'nested'

        # 暂时把frequency < = 0.2的insertion当成 somatic insertion，
        # 后面要根据frequency 的计算方式进行修改
        if self.frequency <= 0.2:
            insertion_type = 'soma|' + nested
        else:
            insertion_type = 'germ|' + nested

        #if self.state:
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.ref_name, self.bp[0]+1, self.bp[-1]+1, self.te_type, 0, strand, self.supp_reads, self.frequency, insertion_type, self.tsd, self.insert_seq_info, self.consensus )
        #else:
        #    continue

    cpdef get_consensus(self, genome_fa,  genome_idx):
        """ define consensus sequence """
        cdef:
            str consensus_prefix  = self.out_path + self.ref_name + "/" + self.c_id
            list insert_size_list=[] # 找到一个合适的组装的seq size
            int i
            
        try:
            consensus_with_tsd = get_consensus(self.seg_list, self.out_path, self.ref_name, self.te_idx, self.c_id, self.type_list, genome_fa,  genome_idx )
        except(IndexError):
            consensus_with_tsd = ["IndexError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
 
        except(KeyError):
            consensus_with_tsd = ["KeyError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        except(ValueError):
            consensus_with_tsd = ["ValueError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        except(UnboundLocalError):
            consensus_with_tsd = ["UnboundLocalError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        except(FileNotFoundError):
            consensus_with_tsd = ["FileNotFoundError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        except(TypeError):
            consensus_with_tsd = ["TypeError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None']
            rm_wtdbg2_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
            rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        self.consensus = consensus_with_tsd[0]
        self.tsd = consensus_with_tsd[1]
        self.insert_seq_info = consensus_with_tsd[2]
        

    def find_tsd(self):
        """ find TSD sequence """
        pass

# MrBleem


## consensus
cdef tuple get_consensus(list seg_list, str out_path, str ref_name, str te_idx, str c_id, list te_type_list, genome_fa,  genome_idx  ):
    cdef:
        int i = 0, l =  100
        list insert_size_list = []
        str consensus_seq_temp_file =  out_path + ref_name + "/" + c_id  + ".consensus.temp.fa"


    consensus_seq_temp_file_ob = open(consensus_seq_temp_file, 'w')

    insert_type_list = list(set([ insert_seq.segname[-1] for insert_seq in seg_list ])) # 为了判断是否包含span reads
    # 在用wtdbg2组装的时候，发现如果组装长度（大致等于插入TE序列长度，这里暂时用insert seq的平均值代替）大于3000的时候，flanksize取200比较好
    # 如果小于3000，flanksize取300比较好
    # 但现在看到的也是个别例子，可能还得细化
    insert_size_list = [ len(insert_seq.seq_t) for insert_seq in seg_list ]
    insert_size = int(sum(insert_size_list) / len(insert_size_list))
    # print(sorted(insert_size_list))


    # 只有clip supporting reads的insertion先跳过，因为效果不好，还在尝试参数
    # 实在不行就跳过
    if len(insert_type_list) == 1 and 'm' not in insert_type_list:
        consensus_seq_print = 'None'
        tsd_print = 'None'
        te_break_info_print = 'None'
        return consensus_seq_print, tsd_print ,te_break_info_print 

    for insert_seq in seg_list :
        insert_type_list = []
        if insert_seq.segname[-1] != 'c':
            i = i + 1
            #consensus_seq_temp_file_ob.write( ">" + str(i) + "\n" + insert_seq.seq_t + "\n" )

            # 因为trim suporting reads是在找 candinate insertion的时候，不好根据整体情况判断falksize该取多少
            # 所有trim的时候flanksize先取300，在这里进行判断是否要再裁减一下
            #if insert_size > 3500:
            # if insert_size <= 0 :
            #    if insert_seq.segname[-1] == 'm':
            #        consensus_seq_temp_file_ob.write( ">" + insert_seq.segname + "\n" + insert_seq.seq_t[l:-l] + "\n" )
            #    if insert_seq.segname[-1] == 'l':
            #        consensus_seq_temp_file_ob.write( ">" + insert_seq.segname + "\n" + insert_seq.seq_t[l: ] + "\n" )
            #    if insert_seq.segname[-1] == 'r':
            #        consensus_seq_temp_file_ob.write( ">" + insert_seq.segname + "\n" + insert_seq.seq_t[:-l ] + "\n" )
            #else:
                
            consensus_seq_temp_file_ob.write( ">" + insert_seq.segname + "\n" + insert_seq.seq_t + "\n" )

    consensus_seq_temp_file_ob.close()

    if max(insert_size_list) >=10000:
        consensus_seq = get_consensus_seq('wtdbg2',consensus_seq_temp_file, out_path, ref_name, c_id )
    else:
        consensus_seq = get_consensus_seq('mafft',consensus_seq_temp_file, out_path, ref_name, c_id ) # ./consensus.py/mafft

    
    # 将consensus seq比对到TE consensus上
    # splice mode可以比对有结构变异的insertion，cigar中表示为N
    # test_inser_bam = align_mm2("/home/boxu/temp/wtdbg2/insert_seq/maff.insert.mmi", out_path + ref_name + "/insert.consensus.fa",  out_path + ref_name, thread_mm2=10, preset="splice")
    #test_inser_bam = align_mm2(te_idx, out_path + ref_name + "/insert.consensus.fa",  out_path + ref_name, thread_mm2=10, preset="splice")

    #test_inser_bam_read = AlignmentFile(test_inser_bam, "rb")
    # g = open("/home/boxu/temp/wtdbg2/for_qc/wtdbg2_divergency.txt", 'a')
    # for read in test_inser_bam_read:
    #    cigar_count = AlignedSegment.get_cigar_stats(read)[0];n_M = cigar_count[0];NM = cigar_count[-1];PM = n_M - NM;
    #    divergence = (n_M - PM ) / n_M
        #g.write(str(divergence)+"\n")
    # g.close()

    consensus_insert_bam = align_mm2(te_idx, out_path + ref_name + "/" + c_id + ".consensus.fa", out_path + ref_name, thread_mm2=10, preset="splice")
    # consensus_insert_bam = align_mm2("/data/tusers/boxu/annotation/dm3/dm3.transposon_for_simulaTE.mmi", out_path + ref_name + "/" + c_id + ".consensus.fa", out_path + ref_name, thread_mm2=10, preset="splice")

    consensus_insert_bam_read = AlignmentFile(consensus_insert_bam, "rb")
    seq_break_list_all = []
    for insert_seq_read in consensus_insert_bam_read:
        cigar_count = AlignedSegment.get_cigar_stats(insert_seq_read)[0];n_M = cigar_count[0];NM = cigar_count[-1];PM = n_M - NM;divergence = (n_M - PM + 0.01 ) / (n_M + 0.01);print(divergence)
        if insert_seq_read.reference_name not in te_type_list:
            continue

        insert_te = insert_seq_read.reference_name

        # print(insert_seq_read.query_name, insert_seq_read.reference_start, insert_seq_read.reference_end, insert_seq_read.reference_length)


        # 通过比对结果两边clip情况判断consensus seq哪一段是TE
        cigar_list = insert_seq_read.cigartuples
        clip_genome_l = cigar_list[0]
        clip_genome_r = cigar_list[-1]
        break_l = 0
        break_r = 0
        if clip_genome_l[0] not in (4,5):
            continue
        elif clip_genome_l[1] < 50:
            continue

        if clip_genome_r[0] not in (4,5):
            continue
        elif clip_genome_r[1] < 50:
            continue

        if clip_genome_l[0] in (4, 5):
            if insert_seq_read.is_reverse:
                break_l = len(consensus_seq) - clip_genome_l[1]
            else:
                break_l = clip_genome_l[1]
        if clip_genome_r[0]  in (4,5):
            if insert_seq_read.is_reverse:
                break_r = clip_genome_r[1]
            else:
                break_r = len(consensus_seq) - clip_genome_r[1]
        
        break_1 = min([break_l,break_r])
        break_2 = max([break_l,break_r])
        
        if insert_seq_read.is_reverse:
            break_seq_strand = '-'
        else:
            break_seq_strand = '+'
        
        # internal deletion

        indel_cigar = insert_seq_read.cigarstring
        if 'N' in indel_cigar : 
            del_pos_list = [int(n_cigar) for n_cigar in re.findall('[0-9]{1,}',indel_cigar.split('N')[0].split('S')[-1])]

            del_len = del_pos_list[-1]
            del_pos = sum(del_pos_list[:-1])
            seq_break_list_all.append((break_1, break_2,":".join([insert_te , str(insert_seq_read.reference_start),str(insert_seq_read.reference_end) ,break_seq_strand,"with_"+str(del_len)+"_del",str(insert_seq_read.reference_start + del_pos), str(insert_seq_read.reference_start + del_pos + del_len)])))

        else:
            seq_break_list_all.append((break_1, break_2,":".join([insert_te , str(insert_seq_read.reference_start),str(insert_seq_read.reference_end) ,break_seq_strand ])))
    seq_break_list_all_sorted = sorted(seq_break_list_all, key = itemgetter(1))

    # print(seq_break_list_all_sorted)


    consensus_seq_print_list = []
    te_break_info_print_list = []
    for i in range(len(seq_break_list_all_sorted)):
        break_p = seq_break_list_all_sorted[i]
        if i == 0:
            consensus_seq_print_list.append(consensus_seq[0:break_p[0]])
            te_break_info_print_list.append("0-"+str(break_p[0])+":genome")
        consensus_seq_print_list.append(consensus_seq[break_p[0]:break_p[1]])
        te_break_info_print_list.append(str(break_p[0])+"-"+str(break_p[1]) + ":" + break_p[2])

        if i == len(seq_break_list_all_sorted) - 1 :
            consensus_seq_print_list.append(consensus_seq[break_p[1]:])
            te_break_info_print_list.append( str(break_p[1]) + "-" + str(len(consensus_seq)) + ":genome" )

    consensus_seq_print = '__'.join(consensus_seq_print_list)
    te_break_info_print = '__'.join(te_break_info_print_list)
    seq_break_list_tsd = [seq_break_list_all_sorted[0][0], seq_break_list_all_sorted[-1][1]]

    
    tsd_print = get_tsd( consensus_seq, seq_break_list_tsd, out_path, ref_name, c_id, genome_fa,  genome_idx )
    #tsd_print = 'None'
    
    # 去掉中间文件，如果保留的话，pysam在读consensus.fa的时候会出问题，感觉是因为内存释放的问题，
    rm_wtdbg2_cmd = ['rm', out_path + ref_name +  "/*" + c_id + "*" ]
    # rm_wtdbg2_proc = Popen([" ".join(rm_wtdbg2_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
    
    
    return  consensus_seq_print , tsd_print , te_break_info_print


cdef str get_tsd(str consensus_seq, list seq_break_list, str out_path, str ref_name, str c_id, str genome_fa,  str genome_idx ):
    tsd_temp_file = out_path + ref_name + "/"+ c_id + ".tsd.temp.fa"
    tsd = ''

    # 取consensus seq两端非TE的序列写入文件与genome进行比对
    tsd_l = consensus_seq[: seq_break_list[0]]
    tsd_r = consensus_seq[seq_break_list[1]: ]
    tsd_temp_file_ob = open(tsd_temp_file, 'w')
    tsd_temp_file_ob.write('>left_tsd\n' + tsd_l + '\n>right_tsd\n' + tsd_r)
    tsd_temp_file_ob.close()
    tsd_bam = align_mm2(genome_idx, tsd_temp_file, out_path + ref_name, thread_mm2=10, preset="map-ont")

    tsd_bam_ob = AlignmentFile(tsd_bam, 'rb')
    tsd_temp_dic = {}
    tsd_temp_list = []
    for read in tsd_bam_ob:
        if read.reference_name != ref_name:
            continue
        cigar_list = read.cigartuples
        tsd_temp_dic[read.query_name] = [read.reference_start, read.reference_end]
        if read.query_name == 'left_tsd':
            if cigar_list[-1][0] not in (4,5):
                left_tsd_clip = 0
            else:
                left_tsd_clip = cigar_list[-1][1]
                # print(left_tsd_clip)
        if read.query_name == 'right_tsd':
            if cigar_list[0][0] not in (4,5):
                right_tsd_clip = 0
            else:
                right_tsd_clip = cigar_list[0][1]
    
    # print(tsd_temp_dic)
    tsd_position = sorted([ tsd_temp_dic['right_tsd'][0],tsd_temp_dic['right_tsd'][1] ,tsd_temp_dic['left_tsd'][0], tsd_temp_dic['left_tsd'][1]])[1:3]
    tsd_len = tsd_position[1] - tsd_position[0]
    # print(tsd_position)
    
    # print(-tsd_len-left_tsd_clip,-left_tsd_clip)
    if tsd_len > 50:
        tsd = "too_long"
        return tsd
    if left_tsd_clip == 0:
        tsd_left = tsd_l[ -tsd_len: ]
    else:
        tsd_left = tsd_l[ -tsd_len-left_tsd_clip: -left_tsd_clip]
    tsd_right = tsd_r[right_tsd_clip : right_tsd_clip + tsd_len]

    tsd_genome = FastaFile(genome_fa).fetch(ref_name, tsd_position[0], tsd_position[1])


    tsd = tsd_left + '__' + tsd_right + '__' + tsd_genome + '__pass'
    # print(tsd)
    return tsd


## trim ##

cdef tuple trim(bint is_reverse, int qlen, str seq, int q_st, int q_en, int op, int op_len, int flanksize):
    """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。 """
    cdef:
        int st_t, en_t
        int q_st_t, q_en_t
        str seq_t

    flanksize_s = flanksize
    flanksize_e = flanksize
    # 如果是clip reads，clip的部分完全保留
    if op in (4, 5):
        if q_st == 0: # l
            flanksize_s = 0
        else:
            flanksize_e = 0

    st_t = q_st - flanksize_s
    en_t = q_en + flanksize_e

    if st_t < 0:
        st_t = 0
    if en_t > qlen:
        en_t = qlen

    q_st_t = q_st - st_t
    q_en_t = q_st_t + op_len
    seq_t = seq[st_t:en_t]


    return seq_t, q_st_t, q_en_t

## extract_seg ##
cdef list extract_seg(object read, object seg_fa, dict read_seq_dic, int flanksize=300, unsigned int min_len=200, unsigned int max_len=11000):
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
        object t


    # MrBleem
    # for hard clip

    # 根据primary alignment序列信息获得suplymentary alignment序列信息
    # read_seq_dic里存储的是read的原始序列信息
    # reads比对到genome的正链，bam文件中的query sequency与raw data fasta文件中的read sequence是一致的
    # 如果是比对到genome的负链，bam文件中的序列与read sequence是反向互补的
    if read.is_supplementary:        
        if read.query_name in read_seq_dic.keys():
            # seq = read_seq_dic[read.query_name]
            if read.is_reverse:
                seq = revcomp(read_seq_dic[read.query_name])
            else:
                seq = read_seq_dic[read.query_name]
            qlen = len(seq)

    ap[-1] = None
    ap[qlen] = None

    # 找insertion通过遍历cigar字段
    # 找到所有可能存在的insertion位置
    # 然后记录下来insert片段和比对到reference的位置
    
    for op, op_len in read.cigartuples:
        if op in (1, 4, 5) and op_len >= min_len: # I or S or H
            if op_len > max_len:
                continue

            if op != 1 and read.query_length < 300:
                continue

            # find query start & end of the segment
            q_st = start_idx    
            q_en = start_idx + op_len
            

            # find reference start & end of the segment
            
            r_st = ap[q_st - 1]

            if op == 5:
                if start_idx == 0:
                    r_st = ap[0]
                r_en = r_st
            else:
                r_en = ap[q_en]
            # if read.query_name == ""


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
            elif op == 5: # H
                if q_st == 0:
                    segname = '{}_{}_{}_{}_{}_{}_H_l'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)
                else:
                    segname = '{}_{}_{}_{}_{}_{}_H_r'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)

            # trim the query sequence and define the trimmed query start & end of the segment
            if seq:
                seq_t, q_st_t, q_en_t = trim(read.is_reverse, qlen, seq, q_st, q_en, op, op_len, flanksize)
                
                # write out trimmed sequence
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, q_st_t, q_en_t, seq_t)
                seg_fa.write('>{}\t{}\n{}\n'.format(segname, len(seq_t), seq_t))
            else:
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, 0, 0)

            segs.append((segname, seg))
            # print(segname)
        
        if op in (0, 1, 4, 7, 8):
            start_idx += op_len

    return segs # list of tuples, (segname, seg)


## align_mm2 ##
cdef str align_mm2(str ref, str query, str outpath, int thread_mm2, str preset="map-ont"):
    cdef:
        int exitcode
        ## str wk_dir = os.path.dirname(query)
        str wk_dir = outpath
        str q_p = os.path.basename(query).rsplit('.', 1)[0]
        str ref_p = os.path.basename(ref).rsplit('.', 1)[0]
        str out_path = os.path.join(wk_dir, ref_p + '_' + q_p + '.bam')
        # list cmd_mm2 = ['minimap2', '-t', str(thread_mm2), '-aYx', preset, '--secondary=no', ref, query, '|', 'samtools', 'view', '-bhSG 4', '-', '|', 'samtools', 'sort', '-o', out_path, '-']
        list cmd_mm2 = ['minimap2', '-t', str(thread_mm2), '-aYx', preset, '--secondary=no',ref, query, '|', 'samtools', 'view', '-bhS', '-', '|', 'samtools', 'sort', '-o', out_path, '-']
        list cmd_idx = ['samtools index ' + out_path]
    
    # pysam.sort("-o", "output.bam", "ex1.bam")

    align_proc = Popen([" ".join(cmd_mm2)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = align_proc.wait()
    if exitcode != 0:
        raise Exception("Error: minimap2 alignment for {} to {} failed".format(q_p, ref_p))
    
    idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = idx_proc.wait()
    if exitcode != 0:
        raise Exception("Error: samtools index for {} failed".format(out_path))
    return out_path


cdef score_te_alignment(object read, dict te_size_dict):
    cdef:
        int score_te_align = 0, map_te_len, insert_te_len
        str map_te_name
    
    
    map_te_name = read.reference_name.split(':')[0]
    map_te_len = int(te_size_dict[map_te_name])
    insert_te_len = read.reference_length
    
    # 计算divergency
    # mismatch / （ mismatch + perfect match ）
    cigar_count = AlignedSegment.get_cigar_stats(read)[0]
    n_M = cigar_count[0]
    NM = cigar_count[-1]
    PM = n_M - NM
    divergence = (n_M - PM ) / n_M


    # 这里考虑read上比对到TE的长度，如果太短就不太可靠，以及如果比对到TE的部分只占read的一小部分，也不太可靠
    mappbility_for_read = float(read.query_alignment_length) / read.query_length
    mappbility_for_te = float(insert_te_len) / map_te_len

    if read.query_name[-1] == 'm':
        
        if mappbility_for_te >= 0.25:
            score_te_align = 1
        #print('ck3\n'+map_te_name,read.query_name)
        #print(read.query_alignment_length,read.query_length,read.reference_start,read.reference_end,read.pos)
        #print(read.query_alignment_start)
        #print(mappbility_for_read, mappbility_for_te)
    else:
        
        mappbility_for_read = float(read.query_alignment_length) / read.query_length
        mappbility_for_te = float(insert_te_len) / map_te_len
        
        #print('ck3\n'+map_te_name,read.query_name)
        #print(read.query_alignment_length,read.query_length,read.reference_start,read.reference_end,read.pos)
        #print(read.query_alignment_start)
        #print(mappbility_for_read, mappbility_for_te)

        if mappbility_for_read > 0.6 and read.mapq >= 40:
            score_te_align = 1
        elif read.pos <= 20:
            if mappbility_for_te >= 0.1 and read.mapq >= 30:
                score_te_align = 1
        else:
            if mappbility_for_te >= 0.25 and read.mapq >= 60:
                score_te_align = 1
            
    
    return score_te_align




## parse_te ##
cdef parse_te_aln(dict seg_dict, str te_bam, dict te_size_dict):
    cdef:
        object te_aln = AlignmentFile(te_bam, "rb"), r, te_not_aln
        str seg_n, te_not_aln_bam = te_bam[:-4] + '.not_aln.bam'
        int q_st, q_en, score
        list uncorect_list = []
    
        list cmd_te_not_aln = ['samtools', 'view', '-bhSf 4', te_bam, '|', 'samtools', 'sort', '-o', te_not_aln_bam, '-']
        list cmd_sam_index = ['samtools index ' + te_not_aln_bam]
    
    # 这里考虑的是比对到reference TE的情况
    # 但后面发现也存在因为insert片段太短而比对不到TE上的，可能还需要调整比对的参数

    ## 把没有比对到TE上的reads取出来
    ## 然后判断这些reads是不是clip的片段
    ## 如果是，那就代表这些clip的片段没有比对到TE上，
    ## 如果真是一个insertion的话，是应该可以比对到TE上的，因为没有clip的部分已经比对到基因组上了
    te_not_aln_proc = Popen([" ".join(cmd_te_not_aln)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = te_not_aln_proc.wait()
    if exitcode != 0:
        raise Exception("Error")
    
    sam_index_proc = Popen(cmd_sam_index, stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = sam_index_proc.wait()
    if exitcode != 0:
        raise Exception("Error: samtools index for {} failed".format(te_not_aln_bam))
    
    te_not_aln = AlignmentFile(te_not_aln_bam, "rb")

    for r in te_not_aln:
        # print(r.query_name[-1])
        if r.query_name[-1] in ('l', 'r'):
            raw_read_id = "_".join(r.query_name.split('_')[3:5])
            uncorect_list.append(raw_read_id)

    # print(uncorect_list)

    for r in te_aln:
        raw_read_id = "_".join(r.query_name.split('_')[3:5])
        if not r.is_unmapped and r.flag !=4 and raw_read_id not in uncorect_list:
            score_te = score_te_alignment(r, te_size_dict)
            if score_te == 1:
                seg_n = r.query_name
                if r.is_reverse:
                    q_st = r.query_length - r.qend
                    q_en = q_st + r.qlen           
                else:
                    q_st = r.qstart
                    q_en = r.qend

                seg_dict[seg_n].add_te(r.reference_name, q_st, q_en, r.reference_start, r.reference_end, int(r.is_reverse))


## collect_seg ##
cdef dict collect_seg(object ref_aln, str chrom, str te_idx, str outpath, dict read_seq_dic, dict te_size_dict, int flanksize, int flank=200):
    cdef:
        object seg_fa = open(outpath+"/"+chrom+".tmp.fa", "w")
        dict seg_dict = {}
        list segs = []
        object read
        str te_bam


    # extract segment from reference genome alignment
    for read in ref_aln.fetch(chrom):
        if not read.is_secondary:
            if read.mapping_quality > 0:
                segs.extend(extract_seg(read, seg_fa, read_seq_dic, flanksize))

    seg_fa.close()
    seg_dict = dict(segs) # should we remove segs after this?

    # align segment sequence to transposon consensus sequence
    te_bam = align_mm2(te_idx, outpath + "/" + chrom + ".tmp.fa", outpath, thread_mm2=10)
    # 这里比对到TE之后，还需要考虑一个insert片段可能会比对到两个以上的的TE consensus上

    # parse TE alignment of segment sequence
    parse_te_aln(seg_dict, te_bam, te_size_dict)

    return seg_dict


## build_cluster ##
cpdef dict build_cluster(str ref_bam, str chrom, str te_idx, str outpath, dict read_seq_dic, dict te_size_dict, int flanksize, int flank=200 ):
    cdef:
        str c_id
        list c_ids
        dict seg_dict
        dict cluster_dict = {}
        object tree = Intersecter(), ref_aln = AlignmentFile(ref_bam, 'rb')
        Segment seg
        Cluster cluster
    
    seg_dict = collect_seg(ref_aln, chrom, te_idx, outpath, read_seq_dic, te_size_dict, flanksize, flank)

    
    for seg in seg_dict.values():
        if seg.n_aln:
            c_ids = [x.value for x in tree.find(seg.rpos-1, seg.rpos+1)]
            if c_ids:
                for c_id in c_ids:
                    cluster_dict[c_id].add_seg(seg)
            else:        
                c_id = str(uuid4())
                cluster = Cluster(c_id, chrom, outpath, ref_aln, te_idx)
                cluster_dict[c_id] = cluster
                cluster_dict[c_id].add_seg(seg)
                tree.add_interval(Interval(seg.rpos-flank, seg.rpos+flank, value=c_id))
        
    return cluster_dict



## process_cluster ##
cpdef process_cluster(dict cluster_dict, str chrom, str out_path, str genome_fa, str genome_idx, str repeatmasker_file):
    cdef:
        str c_id
        Cluster cluster
        object fout = open("{}/{}.tmp.bed".format(out_path,chrom), "w")

    
    mkdir_cmd = ['mkdir', out_path + "/temp/" + chrom ]
    mkdir_cmd_proc = Popen([" ".join(mkdir_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = mkdir_cmd_proc.wait()

    
    for c_id in cluster_dict:
        cluster_dict[c_id].find_bp()
        cluster_dict[c_id].te_stat(repeatmasker_file)
        cluster_dict[c_id].static()
        cluster_dict[c_id].get_consensus(genome_fa, genome_idx)
        # print(cluster_dict[c_id].state)
        if cluster_dict[c_id].state:
            fout.write(cluster_dict[c_id].cluster2bed())

    fout.close()
    return cluster_dict
