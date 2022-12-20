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
#    insert 序列长度在3000以下的，用mafft组装，3000以上的，用wtdbg2组装
#       组装之后的序列可以得到结构信息
#       用raw reads跟TE比对的结果做consensus，得到与insert序列中TE数量对应的consensus序列
#           为了提高准确率，加入reference based 的 polishing / assembly
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
        list type_list, cord_list, orients, score_list
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
        self.score_list = []
        self.read = read
    
    cpdef add_te(self, str te_n, int q_st, int q_en, int r_st, int r_en, int te_strand, float score):
        
        self.cord_list.append( (q_st, q_en, r_st, r_en, te_strand) )
        self.orients.append( te_strand )
        self.score_list.append(score)

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
        str te_idx, divergency, consensus_meth, su_type, te_size, genome_fa
        list score_list


    def __init__(self, str c_id, str ref_name, str out_path, object ref_aln, str te_idx, str genome_fa):
        self.c_id = c_id
        self.ref_name = ref_name
        self.out_path = out_path
        self.seg_list = []
        self.type_list = []
        self.orients = []
        self.score_list = []
        self.supp_reads_list = []
        self.unsupp_reads_list = []
        self.span_reads_list = []
        self.consensus = ''
        self.tsd = ''
        self.insert_seq_info = ''
        self.state = 1 # insertion 的是否保留的标签，现在为了分清每个insertion的具体情况，有多种state，会保留下来的stat有 1、16、10、12、14、24 
        self.ref_aln = ref_aln
        self.genome_fa = genome_fa
        self.frequency = 0
        self.supp_reads_num = 0
        self.segname_list = []
        self.te_idx = te_idx
        self.consensus_meth = 'none'

        self.divergency = '0'    # check
        self.su_type = 'none'   # check  # supporting reads type ，有三类，完全跨过insertion，两端clip，测试用的，后面会去掉
        self.te_size = '0' # check  # 看insert size有多大

    cpdef add_seg(self, Segment seg):
        """ 向Cluster对象的segs当中添加周围的Segment """
        self.seg_list.append(seg)
        self.type_list.extend(seg.type_list)
        self.orients.extend(seg.orients)
        self.supp_reads_list.append(seg.read)
        self.segname_list.append(seg.segname)
        self.score_list.extend(seg.score_list)

    cpdef find_bp(self):
        """ find break points """
        cdef:
            list bps
            Segment s
        
        bps = sorted([s.rpos for s in self.seg_list]) # Q22: left/right breakpoint 是否要分开找呢？因为不同类型的 segment 理论上能提供的 breakpoint 类型也不同
        self.bp = [bps[0], bps[-1]]


        

    cpdef te_stat(self, repeatmasker_file, te_size_dict):
        """
        用于判断该 insertion cluster 该保留还是丢弃，根据 cluster 区间与 Repeat masker 注释区间是否有 overlap
        来进行判断，有则丢弃，无则保留

        Parameters:
            repeatmasker_file: string
            Repeatmasker 注释结果（BED文件）.

        Returns：
            self.state  
                0代表丢弃，1代表保留
      
        """
        cdef:
            str te_temp
            list te_temp_list=[]
        # 暂时输出一个insertion位置上所有可能的TE type
        # 但nested insertion的方向还没处理

        # 写入文件跟rmk对比去掉repeatmasker附近的insertion

        
        # 将一个insertion cluster里找到的可能的insertion结果写入文件
        # TE1   insertion_start    insertion_end
        # TE2   insertion_start    insertion_end
        # ....
        g = open( self.out_path + self.ref_name + "/" + self.c_id + "te_type.bed", 'w') # Q23: 因为 self.out_path 来源于 SMS.py 当中定义的 out_path，所有有可能末尾不带 '/'

        for te_temp in list(set(self.type_list)): # EXP: type_list: [tename1_|_0, tename2_|_1, ...]
            te_temp_list.append(te_temp.split('_|_')[0]) # EXP: te_temp_list: [tename1, tename2, ...]
            self.orient = int(te_temp.split('_|_')[1]) # Q24: 此处直接修改了 cluster 的 orient 属性，后面是否有其他的修改？

            if self.orient == 1:
                te_temp_strand = '-'
            else:
                te_temp_strand = '+'
            g.write("\t".join([self.ref_name, str(self.bp[0]-5), str(self.bp[1]+5), te_temp.split('_|_')[0], str(0), te_temp_strand ]) + '\n')

        
        self.te_type = "|".join(te_temp_list) # Q25: 无效操作？因为下面又给这个属性赋予了新的值
        g.close()

        # 前面写入的信息与repeat masker做intersect，保留非repeat masker的部分
        rm_rmk_cmd = ['bedtools', 'intersect', '-a', repeatmasker_file, '-b', self.out_path + self.ref_name + "/" + self.c_id + "te_type.bed" ,'-wa','-wb', '|','awk', '\'{if($4==$10 && $6==$12){print $4}}\'' ] # Q26: 输出结果的每一列是啥内容？此外，对于每一个 cluster，都会重复这一步 1-vs-N 的 intersection
        rmk_te_list = os.popen(" ".join(rm_rmk_cmd)).readlines() # Q27: 这个 bedtools 的结果只有一行吗？
        if len(rmk_te_list) > 0:
            rmk_te = rmk_te_list[0].strip()
            if rmk_te in te_temp_list:
                te_temp_list.remove(rmk_te)
        
        self.type_list = te_temp_list
        self.te_type = "|".join(self.type_list)

        self.te_size = "|".join([str(te_size_dict[te]) for te in self.type_list]) # check
        # type_list 是这个cluster包含的TE insertion类型数量
        # 如果是0，可能是因为这个地方是一个repeat masker，就丢掉这个地方
        # 这里的state后面不会保留
        # Q28: 所以这里 state 的赋值逻辑是什么？看起来 te_stat 这个方法是对原本的 self.type_list 去重、去除与 repeatmasker result overlap 且 te name 相同的元素？
        if len(self.type_list) == 0:
            self.state = self.state + 1

        if len(self.seg_list) == 1 and len(self.type_list) == 1 :
            if self.seg_list[0].score_list[0] < 1 or self.seg_list[0].read.mapping_quality < 50 :
                self.state = self.state + 3

    # MrBleem
    cpdef static(self):
        """ 筛选supporting reads和unsuporting reads，计算frequency """
        """
        Function:
            用于判断该 insertion cluster 该保留还是丢弃，根据 cluster 区间与 Repeat masker 注释区间是否有 overlap
            来进行判断，有则丢弃，无则保留

        Parameters:
            repeatmasker_file: string
            Repeatmasker 注释结果（BED文件）.

        Returns：
            self.state  
                0代表丢弃，1代表保留
                目前可保留的有 10 12 14 16 24
                为了方便测试，暂时保留，后续再做更改
      
        """
        # 这一部分的思路是
        # 先遍历了一遍一个insertion cluster中的所有segment/supporting reads，遍历的时候做了两个事情
        # 1. 合并reads，因为有些segment可能是来自于同一条supporting reads，所以需要根据raw reads id进行一步合并
        #   
        # 2. 判断double clip 的supporting reads
        #   一条reads，如果两端都clip的话，很难说它在genome上比对是正确的，就像船缺少了锚一样，不稳定
        #   如果存在这样的double clip supporting reads，这个insertion就打上标签“black”

        # 根据合并之后的segment list计算supporting reads数目
        #    例如有两个segment：reads1_l, reads1_r, 这样的话，这两个segment只能算一个，只能留下一条insert seq，但这条reads1其实是可以看作是一条spanning reads的，N_supp的数目+2，
        # 同时统计一下supporting reads的类型，n_m : number of spanning reads, n_l: number of left clip reads, n_r : number of right clip reads
        # 如果一个insertion它没有spanning reads，而且只有一边的clip reads，也打上black的标签


        # 根据insertion的大致位置，抓取周围map到genome上的所有reads
        # 在这一步需要做两件事，一件是筛选出来哪些是un supporting reads，一件是根据结构特征判断这个区域的insertion是否可取，或者说，这个区域是否是reference TE
        # 对于un supporting reads筛选条件有跨过insertion 左右breakpoint 30bp(un_support_reads_spand_size) ,这个cutoff理论上应该与op_len的cutoff保持一致，simplify code的时候发现不一致，回头改了测试一下
        #       而且还希望mapping quality高一点 >50，以及divergency <= 0.1
        # 从iGV里看在genome上有些reads比对情况，有明显的两端截断情况，而且很整齐，之前在组会上有展示过
        #    这种情况多见于reference TE以及NNNNNN的区域
        #    所以根据这个结构特征去做一步筛选
        #    判断方法是统计insertion 区域double clip reads有多少条，这些double clip的reads 可能不是supporting reads，也可能不是un supporting reads，但它的比对情况可以让我们了解这个区域的比对结构
        #       但有的insertion虽然supporting reads包含double clip reads，但两边都有，而且有跨过这个insertion的reads(supporting reads和un supporting reads都算)，这部分是可以留下来的


        # 如图所示：
        # 这里面的reads没有箭头的指double clip reads
        #       -------------------------                          -------------------------    -------------------------
        #       -------------------------                          -------------------------    -------------------------
        #       -------------------------                          -------------------------    -------------------------
        #       -------------------------                          ------------------------------------------------------》（spinning reads）
        #      ⬆️                       ⬆️                                                    ⬆️
        #  【insertion1】        【insertion2】                                          【insertion3】

        # 这里insertion1 和insertion2就会被丢掉去，insertion3就会被保留




        # supporting reads
        cdef:
            dict black_region = {'start':100000000000,'end':0}
            int black_read_num = 0, sup_black_read_num = 0, l_black_read_num =0,  r_black_read_num = 0,black_n_l = 0,black_n_r = 0, double_clip_reads_num = 0
            int span_read_left_clip_len = 0, span_read_right_clip_len = 0
            list support_span_reads_list = []
            str span_tag = 'none'
            list supp_reads_list_id = [], supp_reads_m_list_id = [], supp_reads_c_list_id = [], new_seg_list = [], supp_reads_type = []
            int num_supp = 0, un_support_reads_spand_size = 30
            dict  supp_reads_dict = {}

            # self.state insertion 的是否保留的标签，有两个值，0和1，0表示这个insertion是假的，会丢掉，1就会被留下
            # 现在为了分清每个insertion的具体情况，有多种state，会保留下来的stat有 1、16、10、12、14、24 
        
        
        
        if self.state == 0: # Q29: 没有看到 state 会是 0 的情况？
            return 0
        

        black_tag = 'none'

        for candi_seg in self.seg_list:
            candi_seg_id = '_'.join(candi_seg.segname.split('_')[3:5]) # Q30: qname_qstart，对于每个 segment 而言是 unique 的，并且有一些 qname 本身可能就包含 '_' 符号？

            if candi_seg_id not in supp_reads_dict: # Q30: 这一步应该是要填充 supp_reads_dict，qname_qstart -> segname, 但是最终每一个 value 长度都只为 1，因为 key 都是 unique 的？
                supp_reads_dict[candi_seg_id] = []
                supp_reads_dict[candi_seg_id].append(candi_seg.segname)
            else:
                supp_reads_dict[candi_seg_id].append(candi_seg.segname)
            

            if candi_seg.segname[-1] == 'm': # Q31: 循环结束之后，new_seg_list 应该与 self.seg_list 相同，num_supp 应该仍然为 0，因为 supp_reads_m_list_id、supp_reads_c_list_id 一直为 0 ?
                if candi_seg_id in supp_reads_m_list_id:
                    continue
                elif candi_seg_id in supp_reads_c_list_id:
                    num_supp = num_supp + 1
                    continue
                else:
                    new_seg_list.append(candi_seg)
                
            else:
                if candi_seg_id in supp_reads_c_list_id:
                    num_supp = num_supp + 1
                    continue
                elif candi_seg_id in supp_reads_m_list_id:
                    continue
                else:
                    new_seg_list.append(candi_seg)

            
            # 处理reads比对到genome上是两端截断的情况，这种区域定义为black region
            # 如果两端截断都超过100，而且数量比较多，那么可能这个区域是genome上的reference TE
            # 作出这种处理是因为在IGV上看出这种两端截断的结构很明显
            # chr2L:2,932,687-2,936,323

            supp_read = candi_seg.read

            if supp_read.cigartuples[0][0] in (4,5) and supp_read.cigartuples[-1][0] in (4, 5):
                if supp_read.cigartuples[0][1] > 30 and supp_read.cigartuples[-1][1] > 30:
                    print('ck')
                    print(candi_seg.segname)
                    if candi_seg.segname[-1] == 'l': # Q32: l_black_read_num、r_black_read_num 似乎没用？
                        l_black_read_num = l_black_read_num + 1
                    if candi_seg.segname[-1] == 'r':
                        r_black_read_num = r_black_read_num + 1

                    sup_black_read_num = sup_black_read_num + 1
                    
                    if black_region['start'] > supp_read.reference_start:
                        black_region['start'] = supp_read.reference_start
                    if supp_read.reference_end > black_region['end']:
                        black_region['end'] = supp_read.reference_end
                    if sup_black_read_num >= 1: # cutoff 设为1是之前有测试过别的cutoff，1可能比较合理，有一条double clip reads后面就得再进行一次判断
                        black_tag = "black"
                    
            else:
                support_span_reads_list.append(supp_read) # Q32: 'span' 也包括 clip 的 reads?

        supp_reads_list_id = supp_reads_dict.keys() # EXP: [qname1_qstart1, qname2_qstart2, ...]
        print('>supp_reads_list_id')
        print(supp_reads_list_id)




        
        N_supp_2 = 0
        n_m = 0
        n_l = 0
        n_r = 0
        # 在这里是为了合并一条reads分开比对的情况
        for srd in supp_reads_dict:
            print('?')
            print(supp_reads_dict[srd])
            s_temp_list = list(set(supp_reads_dict[srd])) # Q33: [segname]，长度为1？原因和 Q30 一样
            s_temp_list_type = [ s[-1] for s in s_temp_list ]
            n_m_t = s_temp_list_type.count('m')
            n_l_t = s_temp_list_type.count('l')
            n_r_t = s_temp_list_type.count('r')
            if n_m_t > 0:
                N_supp_2 = N_supp_2 + 2
                n_m = n_m + 1
            else:
                if n_l_t >= 1 and n_r_t >= 1: # Q33: 应该不会有这种情况出现？
                    N_supp_2 = N_supp_2 + 2
                    n_m = n_m + 1
                else:
                    N_supp_2 = N_supp_2 + 1
                    if n_l_t == 0:
                        n_r = n_r + 1
                    else:
                        n_l = n_l + 1
        
        



        self.supp_reads = str(len(supp_reads_list_id)) +"_|_"+ "_|_".join(supp_reads_list_id) # EXP: N_|_qname1_qstart1_|_qname2_qstart2
        self.su_type = "_".join([str(n_m),str(n_l),str(n_r)])
        
        
        

        black_n_l = 0
        black_n_r = 0
        if n_m == 0 and sup_black_read_num == 0:
            if n_l == 0 or n_r == 0:
                black_tag = 'black'
                # 如果只有一端clip的reads，并且没有balck reads，black_region就考虑insertion 周围就行了
                black_region['start'] = self.bp[0]
                black_region['end'] = self.bp[-1]

        # unsupporting reads
        flanking_region_size = 200
        if self.bp[0] - flanking_region_size <= 0 :
            start_bp = 1
        else:
            start_bp = self.bp[0] - flanking_region_size
        
        for span_read in self.ref_aln.fetch(contig=self.ref_name, start=start_bp, stop=self.bp[-1] + flanking_region_size):

            if black_tag == 'black':
                # 如果有一条好的reads跨过去，这个区域的map可能是正确的
                # 指比对质量高，两段没有长的clip
                
                if span_read.cigartuples[0][0] in (4,5):
                    span_read_left_clip_len = span_read.cigartuples[0][1]
                else:
                    span_read_left_clip_len = 0

                if span_read.cigartuples[-1][0] in (4,5):
                    span_read_right_clip_len = span_read.cigartuples[-1][1]
                else:
                    span_read_right_clip_len = 0

                # 如果一条spanning reads两边clip 的长度比较短的话，10bp，也可以认为这个insertion区域有reads跨过
                if span_read_left_clip_len < 10 and span_read_right_clip_len < 10:
                    if black_region['start'] - span_read.reference_start > 100 and span_read.reference_end - black_region['end'] > 100:
                        span_tag = 'span'
                        #print(span_tag)
                        if span_read.query_name in supp_reads_list_id:
                            black_tag = 'none'

                #print(span_read.query_name,span_read_left_clip_len,span_read_right_clip_len)
                if span_read_left_clip_len > 30 and span_read_right_clip_len > 30:
                    # 如果insertion是在double clip之间的，两条以上，极大概率是假的
                    # Q34: 下面的判断条件与要求正好相反？
                    if span_read.reference_end - self.bp[-1] > 25 and self.bp[0] - span_read.reference_start > 25:
                        double_clip_reads_num = double_clip_reads_num + 1
                        if double_clip_reads_num >= 2: 
                            self.state = 0
                            print('ck_clip_sides')
                            return 0
                    #print(self.bp)
                    #print(span_read.reference_start, span_read.reference_end)    
                    black_read_num = black_read_num + 1
                    if abs(span_read.reference_end - self.bp[0]) < 20:
                        black_n_r = black_n_r + 1
                    if abs(span_read.reference_start - self.bp[-1]) < 20:
                        black_n_l = black_n_l + 1


                
            # count unspporting reads
            # fetch得到insertion position区域的所有reads
            # 如果reads跨过insertion position两端各30bp以上，而且不是supporting reads，这条reads认为是unsuporting reads
            if span_read.query_name not in supp_reads_list_id:
                if self.bp[0] - span_read.reference_start >= un_support_reads_spand_size and span_read.reference_end - self.bp[-1] >= un_support_reads_spand_size:  
                    divergency_uns = cal_divergency(span_read, 'genome')

                    if span_read.mapping_quality >= 60 and divergency_uns <= 0.1:
                        #print(span_read.query_name)
                        if span_read.query_name not in self.unsupp_reads_list:
                            self.unsupp_reads_list.append(span_read.query_name)

        print(">black_read_num")
        print(sup_black_read_num, black_read_num,black_n_l,black_n_r)
        print(n_m, black_tag, span_tag)

        # 筛选
        # 这里主要筛选有double clip的区域
        # 首先在第一轮遍历中supporting reads中如果有double clip的reads，这个insertion或者区域打上black的标签
        # 如果没有double clip，supporting reads只有right clip或soft clip的reads，也打上black的标签
        # 这样的区域，如果有横跨这个区域或者说包含整个insertion的supporting reads，这个区域就是安全的，可以被保留下来

        # 在black区域，除了spanning supporting reads外，需要考虑的几点有
        # unsupporting spanning reads ｜ 这个在便利span_reads时已经做了统计，如果有，就打上 span的标签，没有就是none
        # sup_black_read_num ｜ supporting reads中double clip的read的数目，同时也有统计right和left的分布，极端一点，如果全是black supporting read，这个insertion不太可信
        # black_read_num ｜ 这个区域所有double clip reads数目，这是为了区分有的insertion虽然没有black supporting read，但这个区域本身就是好多balck reads，可能得丢掉
        # black_n_l, black_n_r, n_l, n_r 这几个对象也是在前几点的基础上的辅助参考，如果只有单边的reads支持，这个insertion很难站得住脚 

        
        # 这里的条件判断一言难尽
        # Q35: 这一步是为了根据 "(right/left) black reads"、"(right/left) black support reads" 的数量来设置 cluster 的state？不太能看懂下面这些判断的含义
        if n_m ==0 and black_tag == 'black':
            if sup_black_read_num != 0:
                if span_tag == 'span':
                    if black_read_num <= 3 :
                        if black_read_num == sup_black_read_num:
                            if (black_read_num == black_n_l and n_r > 0 ) or (black_read_num == black_n_r  and n_l > 0):
                                self.state = self.state + 15 # 可留
                            else:
                                self.state = self.state + 17 
                        else:
                            self.state = self.state + 19 # 
                    else:
                        self.state = self.state + 21 # 丢弃的
                else:
                    if n_l != 0 and n_r != 0 :
                        self.state = self.state + 23 # 可留
                    else:
                        self.state = self.state + 25
            else:
                if black_read_num > 0:
                    self.state = self.state + 27 # 丢弃的
        

        # 为了看insertion周围是不是包含 NNNNN 区域
        desert_region_size = 500
        if self.bp[0] - desert_region_size <= 0 :
            start_bp = 1
        else:
            start_bp = self.bp[0] - desert_region_size
        span_genome_fa = FastaFile(self.genome_fa).fetch(self.ref_name, start_bp, self.bp[-1] + desert_region_size)

        # 这里是处理suporting reads里左右两边clip reads数量比例的问题
        # 如果有一条reads是跨过insertion的，多少clip reads都无所谓
        # 如果没有跨过insertion的reads，且只有一边有比较多的（>=2）clip reads,这个insertion也会丢掉
        # 如果两边都有，则计算比例，0.1883698 是算了一下，所有只有clip reads的insertion，左右两边clip reads的比例，拟合出来一个正太分布
        # 取的是正太分布的（平均值 - 3 * 标准差  ～ 平均值 + 3 * 标准差 ）

        # 最开始只有根据left clip和right clip之间的比例来筛选，后面才加了double clip的部分
        # 加了double clip的判断之后发现能够分辨出只有一端clip supporting reads的insertion是真是假，但因为想再iGV里验证一下，所以对不同的情况都加了state，方便从结果中筛选


        # if n_m == 0 and len(self.unsupp_reads_list) == 0 :
        if n_m == 0:
            if 'NNNNNNNNNNNNNNNNN' in span_genome_fa:
                l_desert_N = span_genome_fa.count('N')
                if l_desert_N > 100:
                    print('>0')
                    self.state = self.state + 7
                    #return 0

            if n_l == 0 :
                if n_r >= 4:
                    self.state = self.state + 9
                    print('>1')
                    #return 0
            elif n_r == 0 :
                if n_l >= 4:
                    print('>2')
                    self.state = self.state + 11
                    #return 0
            else:
                # 0.1883698是统计了insertion的n_l / n_r 的比例，然后做了个分布，拟合了正太分布模型，取 x-3U x+3U 的区域
                if  abs( 0.5 - ( float(n_l) / sum([n_l,n_r])) ) >  0.5 - 0.1883698: 
                    print('>3')
                    self.state = self.state + 13
                    #return 0

        # frequency

        # 最终选择的frequency的计算方法
        # a = 2; b = 1; c = 2
        # n_supp = a * n_m  + b * ( n_l + n_r )
        # n_unsupp = c * len(self.unsupp_reads_list)
        # frequency2 = float(n_supp) / ( n_supp + n_unsupp )


        self.frequency = float(N_supp_2) / ( N_supp_2 + 2*len(self.unsupp_reads_list) )

        self.seg_list = new_seg_list





    cpdef str cluster2bed(self):
        cdef:
            str strand

        if self.orient:
            strand = '-'
        else:
            strand = '+'
        
        is_nested = ''
        if len(list(set(self.type_list))) > 1:
            is_nested = 'nested' # 前面一步因为跟repeatmasker的注释文件做了intersect，type list中去掉了reference TE，所以，如果一个insertion插入到reference TE中，这里是没有标记nested的

        if len(self.seg_list) < 5 and self.frequency <=  2.0 / ((2*len(self.unsupp_reads_list)-1) + 2):
        # if self.frequency < 0.1428: # 1条clip reads ，3条un supporting reads
            insertion_type = 'soma|' + is_nested
            # self.su_type.split('_')[0] 这个是指supporting reads类型中spanning supporting reads的数量
            if int(self.su_type.split('_')[0]) >= 1 and len(self.seg_list) > 1 :
                insertion_type = 'germ|' + is_nested
        else:
            insertion_type = 'germ|' + is_nested
            



        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.ref_name, self.bp[0]+1, self.bp[-1]+1, self.te_type, self.frequency,  strand, self.supp_reads, self.divergency, self.consensus_meth,insertion_type, self.tsd, self.insert_seq_info, self.consensus, self.su_type, self.te_size, self.state )


    cpdef get_consensus(self, genome_fa,  genome_idx, te_anno_fa):
        """ define consensus sequence """
        """
        Function:
            根据supporting reads组装consensus sequence
            根据flanking sequence，比对到genome，找到TSD的位置及序列信息 
            根据insert seq比对TE的结果对TE部分进行polishing
            
        Parameters:
            seg list
            genome fa
            genome idx
        Returns:
            consensus sequence with tsd
        """
        cdef:
            str consensus_prefix  = self.out_path + self.ref_name + "/" + self.c_id
            list insert_size_list=[] # 找到一个合适的组装的seq size
            int i
            
        try:
            consensus_with_tsd = get_consensus(self.seg_list, self.out_path, self.ref_name, self.te_idx, self.c_id, self.type_list, genome_fa,  genome_idx, te_anno_fa )
            if consensus_with_tsd[-1] != 'NA':
                self.bp = [consensus_with_tsd[-2],consensus_with_tsd[-1]]
        except(IndexError):
            consensus_with_tsd = ["IndexError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None', '0', 'None']
        except(KeyError):
            consensus_with_tsd = ["KeyError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None', '0', 'None']
        except(ValueError):
            consensus_with_tsd = ["ValueError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None', '0', 'None']
        except(UnboundLocalError):
            consensus_with_tsd = ["UnboundLocalError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None', '0', 'None']
        except(FileNotFoundError):
            consensus_with_tsd = ["FileNotFoundError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None','0', 'None']
        except(TypeError):
            consensus_with_tsd = ["TypeError", "__".join([self.ref_name, str(self.bp[0]+1), str(self.bp[-1]+1)]), 'None', '0', 'None']

        # 去掉中间文件，如果保留的话，pysam在读consensus.fa的时候会出问题，感觉是因为内存释放的问题，
        rm_temp_c_id_cmd = ['rm', self.out_path + self.ref_name + "/*" + self.c_id + "*"  ]
        rm_temp_c_id_proc = Popen([" ".join(rm_temp_c_id_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )

        self.consensus = consensus_with_tsd[0]
        self.tsd = consensus_with_tsd[1]
        self.insert_seq_info = consensus_with_tsd[2]
        self.divergency = consensus_with_tsd[3]
        self.consensus_meth = consensus_with_tsd[4]
        

    def find_tsd(self):
        """ find TSD sequence """
        pass

# MrBleem


## consensus
cdef tuple get_consensus(list seg_list, str out_path, str ref_name, str te_idx, str c_id, list te_type_list, str genome_fa,  str genome_idx , str te_anno_fa ):
    """
        Function:
            用supporting reads组装一条consensus sequence，
            并比对到genome和transposon consensus sequence上，
            通过比对结果判断 insert sequence 的结构信息
        Parameters:
            seg_list - 一个insertion的cluster信息
            out_path - temp文件路径
            ref_name - insertion 所在染色体名字
            te_idx - transposon consensus sequence 的minimap2 index
            c_id - 每个insertion的特有id，为了记录中间文件用的，这样在并行的时候不会出错
            te_type_list - insert sequence可能包含的TE类型，已经去除过repeatmasker
            genome_fa - genome reference fasta信息，为了获取tsd区域的序列信息
            genome_idx - genome的minimap2 index

        Returns:
            consensus_seq_print - 最终打印在结果里的序列信息，包含genome序列信息，insert seq的序列信息
            tsd_print - tsd在insert consensus sequence 左右两端的序列信息以及genome上的TSD序列信息
            te_break_info_print - 最终consensus sequence 各段分割点信息
            divergency - consensus sequence 跟transposon reference比得到的divergency
            consensus_meth - 做consensus 的软件类型
    """
    cdef:
        int i = 0, l =  100
        list insert_size_list = []
        str consensus_seq_temp_file =  out_path + ref_name + "/" + c_id  + ".consensus.temp.fa"
        float divergency, divergency_all # ......这一步当中的所有类型的divergency，是为了评估组装的效果怎么样，最后应该不会出现在最终的bed文件中，当然放一个mismatch的divergncy也是可以的

    # return 'none', 'none', 'none'
    consensus_seq_temp_file_ob = open(consensus_seq_temp_file, 'w')

    insert_type_list = list(set([ insert_seq.segname[-1] for insert_seq in seg_list ])) # 为了判断是否包含span reads
    
    insert_size_list = [ len(insert_seq.seq_t) for insert_seq in seg_list ]
    # insert size 判断丢掉了，本来是根据不同的insert size决定insert seq上取多长用于local assembly
    # 因为测试的时候发现wtdbg2对于不同长度组装效果不一样
    # 有的短一点组装效果好些，有些长一点组装效果更好
    # 一般是insertion长度3000以下的话，flanking的长度取长一些比较好，3000以上的话，flanking_size取短一些比较好


    
    # print(insert_type_list)
    # 只有clip supporting reads的insertion先跳过，因为效果不好，还在尝试参数
    # 实在不行就跳过
    if len(insert_type_list) == 1 and 'm' not in insert_type_list: # Q37: 如果这里判断成功，会使得 consensus_seq_temp_file_ob 不能正常关闭
        consensus_seq_print = 'None'
        tsd_print = 'None'
        te_break_info_print = 'None'
        divergency = '0'
        consensus_meth = 'None'
        tsd_positions_st = 'NA'
        tsd_positions_en = 'NA'

        return consensus_seq_print, tsd_print ,te_break_info_print, divergency, consensus_meth, tsd_positions_st, tsd_positions_en

    seq_raw_id_uniq_list = []  # 记录read的raw id，有的reads 可能在map时同时有 right clip和 left clip两种情况，如果重复的话，只取其中一条就够了，不然supporting reads会有重复
    for insert_seq in seg_list :
        insert_type_list = [] # Q38: insert_type_list 在这个循环当中完全无用？
        seq_raw_id = "_".join(insert_seq.segname.split('_')[3:5]) # EXP: qname_qstart
        if seq_raw_id not in seq_raw_id_uniq_list:
            seq_raw_id_uniq_list.append(seq_raw_id)
        else:
            continue # Q38: 阴差阳错下，每个 seq_raw_id 都是 unique 的，无法判断至此，但如果判断至此，就会直接跳过本轮循环
        if insert_seq.segname[-1] != 'c': # Q38: 为什么会是 'c' 呢
            i = i + 1 # Q38: i 的作用是什么？
            consensus_seq_temp_file_ob.write( ">" + insert_seq.segname + "\n" + insert_seq.seq_t + "\n" )

    consensus_seq_temp_file_ob.close()


    # 目前取的flanking size在1500
    if max(insert_size_list) < 0:
        consensus_seq = get_consensus_seq('mafft',consensus_seq_temp_file, out_path, ref_name, c_id, 'local' )
        consensus_meth = 'mafft'
    else:
        consensus_seq = get_consensus_seq('wtdbg2', consensus_seq_temp_file, out_path, ref_name, c_id, '1' )
        consensus_meth = 'wtdbg2'
         # ./consensus.py/mafft



    # reference based polishing
    insert_sequence_bam = align_mm2(te_idx, consensus_seq_temp_file, out_path + ref_name, thread_mm2=10,  mismatch_model="--eqx", preset="splice")
    insert_sequence_bam_obj = AlignmentFile(insert_sequence_bam, 'rb') # EXP: cluster 当中，所有 segment trimmed sequence 与 transposon 的比对结果

    te_anno_fa = '/data/tusers/boxu/annotation/dm3/dm3.transposon_for_simulaTE.fa'
    polish_file_name = out_path + ref_name + "/" + c_id + "_polished.fa"
    
    insert_con_sequence_dict = {}

    polish_file = open(polish_file_name, 'w') # check
    # Q39: 基于 segment trimmed sequences-vs-transposon 的 alignment pileup, 对 te_type_list 当中的每一种 transposon 进行 "polishing"? 没有 coverage 的地方怎么处理？
    for te_reference in te_type_list: # EXP: te_type_list: [tename1, tename2, ...], 这里的 te_type_list 是 set 处理过后的
        insert_con_sequence_dict[te_reference] = ''

        insert_con_sequence = ''
        te_seq = FastaFile(te_anno_fa).fetch(te_reference)
        for pileup_col in insert_sequence_bam_obj.pileup(te_reference):
            # pileup_col 是按顺序取出来每个 reference 对应位置上的比对信息
            # get_query_sequences() 可以获得该位置上的每条 read 的碱基信息
            # Q39: 有点好奇 base_col 长什么样的？因为我的环境下安装的 pysam，貌似没有 get_query_sequences 这个方法？
            base_col = pileup_col.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True)
            anchor_base_list = []
            indel_base_list = []
            for base in base_col:
                anchor_base_list.append(base[0].upper())
                if '+' in base:
                    indel_base_list.append(re.findall( '[ATCG]+', base.split('+')[1].upper() )[0] )

            map_base_list_dict = Counter(anchor_base_list)
            map_base_list = [ (x,map_base_list_dict[x]) for x in map_base_list_dict ]
            map_base_list_sorted = sorted(map_base_list, key = itemgetter(1), reverse=True)
            # insert_con_sequence_dict_for_assem[te_reference] = map_base_list_sorted

            if map_base_list_sorted[0][0] != "*":
                anchor_base = map_base_list_sorted[0][0].split('-')[0]
            
                if len(map_base_list_sorted) >= 2:
                    # 如果一个位置上有两个相同数量的碱基，三个都相同的可能性应该比较低吧，这里暂时没有考虑进去
                    # 如果相同数量的碱基中有reference base ，就取reference的碱基
                    if map_base_list_sorted[0][1] == map_base_list_sorted[1][1]:
                        if te_seq[pileup_col.pos] in [map_base_list_sorted[0][0], map_base_list_sorted[1][0]]: # Q39: 如果都跟 ref te 的 base 不一样呢？
                            anchor_base = te_seq[pileup_col.pos]
                            # print(">"+anchor_base)
                # 当比对结果中有大段 deletion 或者说有splice的话，结果会显示为">" "<"，正向和反向
                if anchor_base == '<' or anchor_base == '>':
                    anchor_base = ''
                insert_con_sequence = insert_con_sequence + anchor_base
            
            # indel after this position
            # indel_base_list
            if len(indel_base_list) == 0:
                continue
            map_base_list_dict = Counter(indel_base_list)
            map_base_list = [(x,map_base_list_dict[x]) for x in map_base_list_dict ]
            map_base_list_sorted = sorted(map_base_list, key = itemgetter(1), reverse=True)
            if map_base_list_sorted[0][0] != "*" and map_base_list_sorted[0][1] >= len(base_col) * 0.5:
                insert_con_sequence = insert_con_sequence + map_base_list_sorted[0][0].split('-')[0]

        insert_con_sequence_dict[te_reference] = insert_con_sequence
        polish_file.write(">"+te_reference+"\n"+insert_con_sequence+"\n") # check
    polish_file.close() # check
    # print(insert_con_sequence_dict)

    # test for polished seq
    # print(te_type_list)
    polish_bam = align_mm2(te_idx, polish_file_name, out_path + ref_name, thread_mm2=10, mismatch_model="--eqx", preset="splice")
    # print(polish_bam)
    polish_bam_obj = AlignmentFile(polish_bam, "rb")
    for polish_seq in polish_bam_obj:
        
        if polish_seq.reference_name not in te_type_list:
            continue
        
        cigar_count = AlignedSegment.get_cigar_stats(polish_seq)[0]
        n_PM = cigar_count[7]
        n_MM = cigar_count[8]
        n_ins = cigar_count[1]
        n_del = cigar_count[2]
        n_M = n_PM + n_MM
        n_all = n_M + n_ins + n_del
        if n_M != 0:
            divergency_polish = "{:.5f}".format( n_MM  / n_M )
            divergency_polish_indel = "{:.5f}".format( ( n_ins + n_del )  / n_all )
            divergency_polish_all = "{:.5f}".format( ( n_MM + n_ins + n_del)  / n_all )
        else:
            divergency_polish = 'none'
            divergency_polish_indel = 'none'
            divergency_polish_all = 'none'

    


    # consensus序列比对到 te 确定结构信息
    mis_mo = "--eqx"
    consensus_insert_bam = align_mm2(te_idx, out_path + ref_name + "/" + c_id + ".consensus.fa", out_path + ref_name, thread_mm2=10, mismatch_model=mis_mo, preset="splice")
    consensus_insert_bam_read = AlignmentFile(consensus_insert_bam, "rb") # EXP: 这个是通过 wtdbg2 生成的 trimmed sequence consensus sequence 与 transposon 比对的结果

    seq_break_list_all = []
    te_insert_strand = {}
    # Q40: 确定 wtdbg2_CSS-vs-TE 比对结果当中每条 alignment 在 CSS 当中的 起始/终止位点、比对到哪个 TE、比对方向、有没有大片段的 deletion，而后按照终止位点进行排序。排序结果储存在 seq_break_list_all_sorted 当中？
    for insert_seq_read in consensus_insert_bam_read:
        
        if insert_seq_read.reference_name not in te_type_list:
            continue
        
        cigar_count = AlignedSegment.get_cigar_stats(insert_seq_read)[0]

        if mis_mo == "--eqx":
            n_PM = cigar_count[7]
            n_MM = cigar_count[8]
            n_ins = cigar_count[1]
            n_del = cigar_count[2]
            n_M = n_PM + n_MM
            n_all = n_M + n_ins + n_del
            if n_M != 0:
                divergency = "{:.5f}".format( n_MM  / n_M )
                divergency_indel = "{:.5f}".format( (n_ins + n_del )  / n_all )
                divergency_all = "{:.5f}".format( (n_ins + n_del + n_M)  / n_all )
            else:
                divergency = 'none'
                divergency_indel = 'none'
                divergency_all = 'none'
        else:
            n_M = cigar_count[0]
            NM = cigar_count[-1]
            PM = n_M - NM
            divergency = (n_M - PM ) / n_M


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
                break_l = len(consensus_seq) - clip_genome_l[1] # Q40: 如果这里的 'l' 和 'r' 是针对 consensus 原本的序列而言的话，那么此处应当是 break_r ？
            else:
                break_l = clip_genome_l[1]
        if clip_genome_r[0]  in (4,5):
            if insert_seq_read.is_reverse:
                break_r = clip_genome_r[1]
            else:
                break_r = len(consensus_seq) - clip_genome_r[1] # Q40: 同上，并且计算方式应该是 break_l = clip_genome_r[1] ?
        
        break_1 = min([break_l,break_r])
        break_2 = max([break_l,break_r])
        
        if insert_seq_read.is_reverse:
            break_seq_strand = '-'
        else:
            break_seq_strand = '+'
        te_insert_strand[insert_seq_read.reference_name] = break_seq_strand
        
        # internal deletion

        indel_cigar = insert_seq_read.cigarstring
        if 'N' in indel_cigar : # Q41: CIGAR string 当中的 N 是指什么？
            del_pos_list = [int(n_cigar) for n_cigar in re.findall('[0-9]{1,}',indel_cigar.split('N')[0].split('S')[-1])]

            del_len = del_pos_list[-1]
            del_pos = sum(del_pos_list[:-1]) # Q41: 这里应该是要计算其它 operation 的长度之和？但是从代码上看，这个 '和' 不包括 'S', 但是包括 'D' 与 'H'(如果有的话)？
            seq_break_list_all.append((break_1, break_2,":".join([insert_te , str(insert_seq_read.reference_start),str(insert_seq_read.reference_end) ,break_seq_strand,"with_"+str(del_len)+"_del",str(insert_seq_read.reference_start + del_pos), str(insert_seq_read.reference_start + del_pos + del_len)])))

        else:
            seq_break_list_all.append((break_1, break_2,":".join([insert_te , str(insert_seq_read.reference_start),str(insert_seq_read.reference_end) ,break_seq_strand ])))
    
    seq_break_list_all_sorted = sorted(seq_break_list_all, key = itemgetter(1))
    # print(seq_break_list_all_sorted)


    polished_sequence = insert_con_sequence 
    polished_sequence = []


    consensus_seq_print_list = []
    te_break_info_print_list = []
    # Q42: 将 consensus sequence 片段与 polished TEs sequence 片段进行拼接
    for i in range(len(seq_break_list_all_sorted)):
        break_p = seq_break_list_all_sorted[i]
        break_te_info = break_p[2].split(':')
        break_te_ref = break_te_info[0]
        break_te_start = break_te_info[1]
        break_te_end = break_te_info[2]

        # left genome
        if i == 0:
            consensus_seq_print_list.append(consensus_seq[0:break_p[0]])
            te_break_info_print_list.append("0-"+str(break_p[0])+":genome")
            polished_sequence.append(consensus_seq[0:break_p[0]])

        # insert sequence
        # Q42: 这里缺 elif 和 else 判断，所以第一个 & 最后一个片段会被 append 两次
        consensus_seq_print_list.append(consensus_seq[break_p[0]:break_p[1]])
        te_break_info_print_list.append(str(break_p[0])+"-"+str(break_p[1]) + ":" + break_p[2])
        if te_insert_strand[break_te_ref] == '-': # Q42: 所以 polishing 的思想是将 polished-TE sequence 片段替换原本的 consensus sequence 片段？
            polished_sequence.append(revcomp(insert_con_sequence_dict[break_te_ref])[int(break_te_start):int(break_te_end)])
        else:
            polished_sequence.append(insert_con_sequence_dict[break_te_ref][int(break_te_start):int(break_te_end)])

        # right genome
        if i == len(seq_break_list_all_sorted) - 1 :
            consensus_seq_print_list.append(consensus_seq[break_p[1]:])
            te_break_info_print_list.append( str(break_p[1]) + "-" + str(len(consensus_seq)) + ":genome" )
            polished_sequence.append(consensus_seq[break_p[1]:])

    consensus_seq_print = '__'.join(consensus_seq_print_list) + "......" + "__".join(polished_sequence)
    te_break_info_print = '__'.join(te_break_info_print_list)
    seq_break_list_tsd = [seq_break_list_all_sorted[0][0], seq_break_list_all_sorted[-1][1]]

    
    tsd_info = get_tsd( consensus_seq, seq_break_list_tsd, out_path, ref_name, c_id, genome_fa,  genome_idx )
    tsd_print = tsd_info[0]
    tsd_positions_st = tsd_info[1]
    tsd_positions_en = tsd_info[2]
    #tsd_print = 'None'
    
    
    
    
    # return  consensus_seq_print , tsd_print , te_break_info_print, str(n_MM) + "_" + str(divergency) + "_" + str(divergency_polish) , consensus_meth
    return  consensus_seq_print , tsd_print , te_break_info_print, str(n_MM) + "_" + str(divergency) + "|" + str(divergency_indel) + "|" + str(divergency_all) + "_" + str(divergency_polish)+ "|" + str(divergency_polish_indel) + "|" + str(divergency_polish_all) , consensus_meth, tsd_positions_st, tsd_positions_en


cdef tuple get_tsd(str consensus_seq, list seq_break_list, str out_path, str ref_name, str c_id, str genome_fa,  str genome_idx ):
    """
        Function:
            由consensus序列获取两端非TE的序列，比对到genome上，确定TSD信息
            用在get_consensus()这个函数里
        Parameters:
            consensus_seq - 由软件（wtdbg2或mafft）获得的consensus 序列
            seq_break_list - 最终consensus sequence 各段分割点信息，用于获取两端flanking region
            剩下同上
        Returns:
            tsd - 左右两端以及genome上的tsd序列信息
                tsd_left + '__' + tsd_right + '__' + tsd_genome + '__pass'

    """
    cdef:
        str tsd = ''
        str tsd_l, tsd_r # tsd left and tsd right
        dict tsd_temp_dic = {}
        list tsd_temp_list = [], tsd_position
        int tsd_len
        str tsd_temp_file = out_path + ref_name + "/"+ c_id + ".tsd.temp.fa"
    

    # 取consensus seq两端非TE的序列写入文件与genome进行比对
    tsd_l = consensus_seq[: seq_break_list[0]]
    tsd_r = consensus_seq[seq_break_list[1]: ]

    tsd_temp_file_ob = open(tsd_temp_file, 'w')
    tsd_temp_file_ob.write('>left_tsd\n' + tsd_l + '\n>right_tsd\n' + tsd_r)
    tsd_temp_file_ob.close()
    
    tsd_bam = align_mm2(genome_idx, tsd_temp_file, out_path + ref_name, thread_mm2=10, mismatch_model="", preset="splice")

    tsd_bam_ob = AlignmentFile(tsd_bam, 'rb')

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
    # Q43: 这里是怎么寻找 TSD 的？没太看懂
    tsd_position = sorted([ tsd_temp_dic['right_tsd'][0],tsd_temp_dic['right_tsd'][1] ,tsd_temp_dic['left_tsd'][0], tsd_temp_dic['left_tsd'][1]])[1:3]
    tsd_len = tsd_position[1] - tsd_position[0]
    # print(tsd_position)
    
    # print(-tsd_len-left_tsd_clip,-left_tsd_clip)
    if tsd_len > 30:
        tsd = "too_long"
        return tsd,'NA','NA'
    if left_tsd_clip == 0:
        tsd_left = tsd_l[ -tsd_len: ]
    else:
        tsd_left = tsd_l[ -tsd_len-left_tsd_clip: -left_tsd_clip]
    tsd_right = tsd_r[right_tsd_clip : right_tsd_clip + tsd_len]

    tsd_genome = FastaFile(genome_fa).fetch(ref_name, tsd_position[0], tsd_position[1])


    tsd = tsd_left + '__' + tsd_right + '__' + tsd_genome + '__pass'
    # print(tsd)
    return tsd,tsd_position[0],tsd_position[1]


## trim ##

cdef tuple trim(bint is_reverse, int qlen, str seq, int q_st, int q_en, int op, int op_len, int flanksize):
    """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。 """
    cdef:
        int st_t, en_t
        int q_st_t, q_en_t
        str seq_t

    flanksize_s = flanksize # Q14: 无需单独区分 flanksize_s、flanksize_e
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

    q_st_t = q_st - st_t # EXP: "query start" of the I/S fragment in the trimmed sequence
    q_en_t = q_st_t + op_len # EXP: "query end" of the I/S fragment in the trimmed sequence
    seq_t = seq[st_t:en_t]


    return seq_t, q_st_t, q_en_t

## extract_seg ##
cdef list extract_seg(object read, object seg_fa, dict read_seq_dic, int flanksize=300,  unsigned int min_len=50, unsigned int max_len=11000 ):
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

    # 根据primary alignment序列信息获得suplymentary alignment序列信息
    # read_seq_dic里存储的是read的原始序列信息
    # reads比对到genome的正链，bam文件中的query sequency与raw data fasta文件中的read sequence是一致的
    # 如果是比对到genome的负链，bam文件中的序列与read sequence是反向互补的
    if read.is_supplementary: # Q9: read_seq_dic 来源于 SMS.py 中定义的 read_seq_dic, 目前是一个空的字典, 此判断目前是无用的？
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
            if op_len > max_len: # Q10: 如果是 op_len > max_len 条件满足，就会跳过此轮循环，进行下一轮，导致本轮循环最后面的语句 start_idx += op_len 无法执行，导致 start_idx 计算错误
                continue

            if op != 1 and read.query_length < 300: # Q11: 和 Q10 一样
                continue

            # find query start & end of the segment
            q_st = start_idx
            q_en = start_idx + op_len # Q12: 如果是 left-hard-clip, q_st 应当 小于 0, 而 q_en 应等于 0
            

            # find reference start & end of the segment
            r_st = ap[q_st - 1]

            if op == 5:
                if start_idx == 0: # Q13: 此处，若是 left-hard-clip，则 r_st = r_en = ap[0]? 
                    r_st = ap[0]
                r_en = r_st
            else:
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
            elif op == 5: # H
                if q_st == 0:
                    segname = '{}_{}_{}_{}_{}_{}_H_l'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)
                else:
                    segname = '{}_{}_{}_{}_{}_{}_H_r'.format(read.reference_name, rpos, int(read.is_reverse), read.query_name, q_st, q_en)

            # trim the query sequence and define the trimmed query start & end of the segment
            if seq:
                seq_t, q_st_t, q_en_t = trim(read.is_reverse, qlen, seq, q_st, q_en, op, op_len, flanksize)
                
              
                
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, q_st_t, q_en_t, seq_t)
                if op in (4,5):
                    if q_st == 0:
                        seg_fa.write('>{}\t{}\n{}\n'.format(segname, len(seq_t[:op_len]), seq_t[:op_len])) # 之前insert seq命名中长度不一致，修改一致
                    else:
                        seg_fa.write('>{}\t{}\n{}\n'.format(segname, len(seq_t[-op_len:]), seq_t[-op_len:]))
                else:
                    seg_fa.write('>{}\t{}\n{}\n'.format(segname, len(seq[q_st:q_en]), seq[q_st:q_en]))
                    

            else:
                # 这一步理论上是走不到的
                print('?')
                seg = Segment(read, segname, read.is_reverse, q_st, q_en, r_st, r_en, rpos, 0, 0)

            segs.append((segname, seg))
            # print(segname)
        
        if op in (0, 1, 4, 7, 8):
            start_idx += op_len

    return segs # list of tuples, (segname, seg)


## align_mm2 ##
cdef str align_mm2(str ref, str query, str outpath, int thread_mm2, str mismatch_model="", str preset="map-ont"):
    # mismatch_model
    # 这个参数对genome是“”
    # 对TE是--eqx 

    # preset
    # 对genome是map ont map-pb
    # 对TE是splice
    cdef:
        int exitcode
        ## str wk_dir = os.path.dirname(query)
        str wk_dir = outpath
        str q_p = os.path.basename(query).rsplit('.', 1)[0]
        str ref_p = os.path.basename(ref).rsplit('.', 1)[0]
        str out_path = os.path.join(wk_dir, ref_p + '_' + q_p + '.bam')
        list cmd_mm2 = ['minimap2', '-t', str(thread_mm2), '-aYx', preset, mismatch_model, '--secondary=no',ref, query, '|', 'samtools', 'view', '-bhS', '-', '|', 'samtools', 'sort', '-o', out_path, '-']
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


cdef cal_divergency(object read, str ref):
    """
        Function:
            根据不同的比对模式计算divergency
        Parameter: 
            reads
            reference
        Returns:
            divergency
    """
    cdef:
        int n_PM
        int n_MM
        float divergency
    cigar_count = AlignedSegment.get_cigar_stats(read)[0]

    if ref == 'te':
        n_PM = cigar_count[7]
        n_MM = cigar_count[8]
        n_ins = cigar_count[1]
        n_del = cigar_count[2]
        n_M = n_PM + n_MM
        n_all = n_M + n_ins + n_del

        divergency = float(n_MM) / ( n_PM + n_MM )

    if ref == 'genome':
        n_M = cigar_count[0]
        NM = cigar_count[-1]
        PM = n_M - NM
        divergency = float((n_M - PM )) / n_M

    return divergency


cdef score_te_alignment(object read, dict te_size_dict):
    """
        Function:
            判断一条supporitng reads跟transposon 比对情况，根据几个参数进行打分0/1
        Parameters:
            read - 从 candidate supporting reads中取出来的 insert sequence 比对的TE上的情况
            te_size_dict - 根据 transposon.size 文件建立的记录各 transposon 大小的字典
        Returns:
            score_te_align - 该比对结果的打分判断结果 0 / 0.5 / 1
                0丢弃 0.5、1暂时保留，如果只有一条supporting reads，且score是0.5，则丢弃

    """
    cdef:
        float score_te_align = 0
        int map_te_len, insert_te_len
        str map_te_name
    
    # Q16: 这里加 split 是因为测试的时候使用的 TE CSS 原本是用来跑 TLDR 的? 据他们的描述，他们的软件能够区分差异细微的 sub-family TE insertion
    map_te_name = read.reference_name.split(':')[0] # 比对的TE的类型
    map_te_len = int(te_size_dict[map_te_name]) # consensus TE的长度
    insert_te_len = read.reference_length # 比对到transposon的长度
    divergency = cal_divergency(read, 'te') # insert seq比对到TE上的divergency，但是这个参数没有考虑到筛选条件中，因为插入进去的TE片段是存在比较高divergency的水平的可能的



    # 这里考虑read上比对到TE的长度，如果太短就不太可靠，以及如果比对到TE的部分只占read的一小部分，也不太可靠
    # Q17: 这里的 'read' 本质上是 trimmed sequence，其中可能同时包含很长的 non-TE 片段，以及相对短的 TE 片段，只看比例对较短的TE不公平。所以我觉得可以考虑对 read.query_alignment_length 设置一个 cutoff ?
    insert_query_length = float( read.query_name.split('_')[-2] )- float(read.query_name.split('_')[-3] )
    mappbility_for_read = float(read.query_alignment_length) / insert_query_length # Q18: 分子的计算方式是否要考虑 flanksize, 因为针对 insert(I) fragment 写出的片段是包括两端 flanksize 片段的
    mappbility_for_te = float(insert_te_len) / map_te_len # Q19: 与 Q18 类似

    if read.query_name[-1] == 'm':


        # test.soma.9 之前这个要求 阈值是0.25
        if mappbility_for_te >= 0.15:  # 如果supporitng reads里有跨过这个insertion的reads，而且reads中transposon的序列与consensus序列比例大于0.25，就保留这条alignment
            score_te_align = 1
    
    else:

        if read.mapping_quality < 30:
            if divergency < 0.02:  
                score_te_align = 1

        if mappbility_for_read > 0.6 and read.mapping_quality >= 40:    # 如果reads大部分序列都是属于consensus TE的，这条比对可能是真的
            score_te_align = 1
        elif read.pos <= 20: # Q20: 对于一些 5' truncated insertion (比如 LINE1) 不太公平?
            if mappbility_for_te >= 0.1 and read.mapping_quality >= 30: # 如果与TE比对是从头开始比对，对reads中TE的占比要求可以低一点，如果不是就要求高一点
                score_te_align = 1
        else:
            if mappbility_for_te >= 0.25 and read.mapping_quality >= 60:
                score_te_align = 1
        
        if read.query_alignment_length / insert_query_length < 0.75:
            score_te_align = 0.5
        
    return score_te_align




## parse_te ##
cdef parse_te_aln(dict seg_dict, str te_bam, dict te_size_dict):
    cdef:
        object te_aln = AlignmentFile(te_bam, "rb"), r
        str seg_n
        int q_st, q_en, score
    
    
    # 这里考虑的是比对到reference TE的情况
    # 但后面发现也存在因为insert片段太短而比对不到TE上的，可能还需要调整比对的参数

    ## 把没有比对到TE上的reads取出来
    ## 然后判断这些reads是不是clip的片段
    ## 如果是，那就代表这些clip的片段没有比对到TE上，
    ## 如果真是一个insertion的话，是应该可以比对到TE上的，因为没有clip的部分已经比对到基因组上了

    # simplify 说明
    # 把not align的部分删掉了，clip的部分可能有多个比对情况，如果因为有一种是没有比对上去就丢掉这个clip 片段，可能会丢掉很多真正的insertion片段
    # 而之前考虑到reads比对到reference TE 造成的clip情况，在后面的判断条件中加了进去，在black region部分


    for r in te_aln:
        raw_read_id = "_".join(r.query_name.split('_')[3:5]) # EXP: queryname_qst of each written segment, e.g. read1_100
        if not r.is_unmapped and r.flag !=4: # Q15: 判断条件重复了？
            score_te = score_te_alignment(r, te_size_dict)
            if score_te > 0:
                seg_n = r.query_name
                if r.is_reverse:
                    q_st = r.query_length - r.qend
                    q_en = q_st + r.qlen
                else:
                    q_st = r.qstart
                    q_en = r.qend

                seg_dict[seg_n].add_te(r.reference_name, q_st, q_en, r.reference_start, r.reference_end, int(r.is_reverse), score_te)


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
                temp_segs = extract_seg(read, seg_fa, read_seq_dic, flanksize)
                segs.extend(temp_segs[0])

    seg_fa.close()
    seg_dict = dict(segs) # EXP: segname -> Segment, btw 转换完之后也许应该删除 segs 这个列表？

    # align segment sequence to transposon consensus sequence
    te_bam = align_mm2(te_idx, outpath + "/" + chrom + ".tmp.fa", outpath, thread_mm2=10, mismatch_model="--eqx")
    # 这里比对到TE之后，还需要考虑一个insert片段可能会比对到两个以上的的TE consensus上

    # parse TE alignment of segment sequence
    parse_te_aln(seg_dict, te_bam, te_size_dict)

    return seg_dict


## build_cluster ##
cpdef dict build_cluster(str ref_bam, str genome_fa, str chrom, str te_idx, str outpath, dict read_seq_dic, dict te_size_dict, int flanksize, int flank=200 ):
    # 这里
    # flanksize 是为了后面提取 insert sequence 两边flank region用的
    # flank 是insertion周围多长合并在一起用的
    cdef:
        str c_id
        list c_ids
        dict seg_dict
        dict cluster_dict = {}
        object tree = Intersecter(), ref_aln = AlignmentFile(ref_bam, 'rb')
        Segment seg
        Cluster cluster
    
    # Q6: 在这个地方，如果 ref_aln 是每个 chromosome 上的 alignment, 会更好吗？
    seg_dict = collect_seg(ref_aln, chrom, te_idx, outpath, read_seq_dic, te_size_dict, flanksize, flank)

    
    for seg in seg_dict.values():
        if seg.n_aln:
            c_ids = [x.value for x in tree.find(seg.rpos-1, seg.rpos+1)]
            if c_ids:
                for c_id in c_ids:
                # Q7: 在添加到 cluster 当中之前，是否要对每个 segment 进行筛选？
                    cluster_dict[c_id].add_seg(seg)
            else:        
                c_id = str(uuid4())
                cluster = Cluster(c_id, chrom, outpath, ref_aln, te_idx, genome_fa)
                cluster_dict[c_id] = cluster
                cluster_dict[c_id].add_seg(seg)
                # Q8: flank 该如何取呢？
                tree.add_interval(Interval(seg.rpos-flank, seg.rpos+flank, value=c_id))
        
    return cluster_dict



## process_cluster ##
cpdef process_cluster(dict cluster_dict, str chrom, str out_path, str genome_fa, str genome_idx, str te_anno_fa,  str repeatmasker_file, dict te_size_dict):
    cdef:
        str c_id
        Cluster cluster
        object fout = open("{}/{}.tmp.bed".format(out_path,chrom), "w")

    
    mkdir_cmd = ['mkdir', out_path + "/temp/" + chrom ]
    mkdir_cmd_proc = Popen([" ".join(mkdir_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = mkdir_cmd_proc.wait()

    
    for c_id in cluster_dict:
        cluster_dict[c_id].find_bp()
        cluster_dict[c_id].te_stat(repeatmasker_file, te_size_dict)
        cluster_dict[c_id].static()
        
        if cluster_dict[c_id].state == 0:
            continue 
        cluster_dict[c_id].get_consensus(genome_fa, genome_idx, te_anno_fa)
        fout.write(cluster_dict[c_id].cluster2bed())

    fout.close()
    return cluster_dict
