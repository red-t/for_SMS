import os
import subprocess
from operator import itemgetter
from collections import Counter
from subprocess import Popen, PIPE, STDOUT, DEVNULL
from tkinter import ANCHOR
from pysam import  FastaFile




###          ###
###  wtdbg2  ###
###          ###
def get_consensus_seq_wtdbg2(consensus_seq_temp_file,out_path,ref_name, c_id, model):
    """
        Function:
            用wtdbg2做 local assembly
            
        Parameters:
            consensus_seq_temp_file - 记录一个insertion supporting reads 序列信息的文件
            out_path - 输出结果文件路径
            ref_name - insertion对应reference chrom名称
            c_id - 每个insertion 的特定id
            model - wtpoa-cns做consensus的两种model
        Returns:
            一条consensus序列
    """
    # wtdbg2 - left
    # wtdbg2_contig_cmd = [  'wtdbg2', '-x', 'preset2', '-q', '-p', '1','-k','6','-S','2', '-L', '0', '-g', '0.5k', '-t', '10', '-i', consensus_seq_temp_file, '-fo', out_path + ref_name + ".consensus.temp2" ]        
    # wtdbg2_contig_cmd = [  'wtdbg2 -x preset2 -L 0 -q -p 0 -k 15 -g 3m -S 2 -e 2 -l 256 --ctg-min-length 1000 --node-len 512 --ctg-min-nodes 2', '-t', '10', '-i', consensus_seq_temp_file, '-fo', out_path + ref_name + ".consensus.temp2" ]        
    
    # -L 对参与组装reads长度的限制，第一步筛选
    # -q quiet
    # -p 和 -k 这两个值指kmer，但不是很理解，这里是用推荐的值，测试后感觉效果还行
    # -e 做有向图时要求contig的最低coverage
    # -l alignment长度
    # --ctg-min-length 1000 组装出来contig的长度
    # --node-len 512  有向图每个node的长度，这里至少两个bin
    # --ctg-min-nodes 2    最终contig只有是由几个node组成的
    wtdbg2_contig_cmd = [  'wtdbg2 -x preset2 -L 0 -q -p 0 -k 15 -g 3m -S 2 -e 1 -l 256 --ctg-min-length 1000 --node-len 512 --ctg-min-nodes 2', '-t', '10', '-i', consensus_seq_temp_file, '-fo', out_path + ref_name + "/" + c_id + ".wtdbg2.consensus.temp2" ]        
    wtdbg2_contig_proc = Popen([" ".join(wtdbg2_contig_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = wtdbg2_contig_proc.wait()

    # print(" ".join(wtdbg2_contig_cmd))
    
    # wtdbg2_consensus_cmd = [ 'wtpoa-cns', '-q','-M','2','-X', '-8','-I', '-2','-D', '-4','-c',model,'-C','3','-F','0.6',  '-t', '10', '-i', out_path + ref_name + "/"+c_id +".wtdbg2.consensus.temp2.ctg.lay.gz", '-fo', out_path + ref_name + "/"+c_id +".consensus.fa" ]
    # wtdbg2_consensus_cmd = [ 'wtpoa-cns', '-q', '-M','2', '-X','-8', '-j', '500','-I','-3', '-D','-5', '-H','-5', '-W','150', '-S','1', '-c',model, '-C','3', '-F','0.2', '-t','10', '-i', out_path + ref_name + "/"+c_id +".wtdbg2.consensus.temp2.ctg.lay.gz", '-fo', out_path + ref_name + "/"+c_id +".consensus.fa" ]
    # 目前最好
    # SMS.test
    # wtdbg2_consensus_cmd = [ 'wtpoa-cns', '-q', '-M','3', '-X','-8', '-j', '1000','-I','-3', '-D','-4', '-H','-3', '-W','150', '-S','1', '-c',model, '-A', '-C','3', '-F','0.5', '-N','20', '-t','10', '-i', out_path + ref_name + "/"+c_id +".wtdbg2.consensus.temp2.ctg.lay.gz", '-fo', out_path + ref_name + "/"+c_id +".consensus.fa" ]
    # SMS
    wtdbg2_consensus_cmd = [ 'wtpoa-cns', '-q', '-M','3', '-X','-10', '-j', '1000','-I','-3', '-D','-4', '-H','-3', '-W','150', '-S','1', '-c',model, '-A', '-C','3', '-F','0.5', '-N','30', '-t','10', '-i', out_path + ref_name + "/"+c_id +".wtdbg2.consensus.temp2.ctg.lay.gz", '-fo', out_path + ref_name + "/"+c_id +".consensus.fa" ]

    wtdbg2_consensus_proc = Popen([" ".join(wtdbg2_consensus_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
    
    exitcode = wtdbg2_consensus_proc.wait()


    # print(len(open(out_path + ref_name + ".consensus.temp.wtdbg2_consensus.fa", 'r').readlines()))
    if len(open(out_path + ref_name + "/"+c_id +".consensus.fa", 'r').readlines()) == 0:
        consensus_seq = ''
    else:
        # fa_index
        # fa_index_cmd = ['samtools', 'faidx', out_path + ref_name + ".consensus.temp.wtdbg2_consensus.fa" ]
        # fa_index_proc = Popen([" ".join(fa_index_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )

        consensus_seq_fa = FastaFile(out_path + ref_name + "/"+c_id +".consensus.fa")
        consensus_seq = consensus_seq_fa.fetch('ctg1')
    return consensus_seq



###           ###
###   mafft   ###
###           ###
def mafft(align_seq_file, model):
    ''' use MAFFT to create MSA '''

    out_file_msa = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
    
    # 根据序列的长度进行参数的选择
    #     if
    # args = [ 'mafft', '--randomseed', '1', align_seq_file ]
    # args = [ 'mafft', '--retree', '1','--op', '0', align_seq_file ]
    # args = [ 'mafft', '--retree', '2', align_seq_file ]
    #args = [ 'mafft', '--localpair', '--maxiterate', '1000', align_seq_file ]
    # args = [ 'mafft', '--localpair', '--maxiterate', '1000', '--op', '0',align_seq_file ]
    # args = [ 'mafft', '--genafpair', '--maxiterate','1000', '--op', '1', align_seq_file ]

    # args = [ 'mafft', '--globalpair','--maxiterate', '500', '--op', '0',align_seq_file ]

    if model == "retree":
        args = [ 'mafft', '--retree', '2', '--op', '0.1',align_seq_file ]
    elif model == "local":
        args = [ 'mafft', '--localpair', '--maxiterate', '1000', '--op', '0.1',align_seq_file ]
    elif model == "gene":
        args = [ 'mafft', '--genafpair', '--maxiterate','1000', '--op', '0.1', align_seq_file ]
    elif model == "global":
        args = [ 'mafft', '--globalpair','--maxiterate', '500', '--op', '0.1',align_seq_file ]

    # args = [ 'mafft', '--retree', '2', '--op', '0',align_seq_file ]
 
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = FNULL)

    with open(out_file_msa, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            # if line[0] != '>':
            out_fa.write(line)
    return out_file_msa



def get_seq_ij_count(seq_info, seq_index_list, anchor_point, i ):
    seq_temp_dict = {}
    for k in seq_index_list:
        seq = seq_info[k]
        seq_temp = ''
                        
        for j in range(anchor_point, i+1):
            if seq[j] != '-':
                seq_temp = seq_temp + seq[j]
        if seq_temp not in seq_temp_dict:
            seq_temp_dict[seq_temp] = [1,k]
        else:
            seq_temp_dict[seq_temp][0] = seq_temp_dict[seq_temp][0]+1
            seq_temp_dict[seq_temp].append(k)


    
    seq_ij_max = sorted([(m,seq_temp_dict[m][0],seq_temp_dict[m][1:],len(m)) for m in seq_temp_dict], key = lambda x:(x[1],x[3]), reverse=True)


    if seq_ij_max[0][0] == '-' or seq_ij_max[0][0] == '':
        seq_ij_max = seq_ij_max[1:]
    
    return seq_ij_max




def consensus_seq(sequences, fre_cuoff, type):
    ids = '*'
    if list(sequences.keys()) != [None]:
        ids = ';'.join(list(sequences.keys()))
    
    seq_info = list(sequences.values())
    seq_len = len(seq_info[0])
    consensus_seq = []
    consensus_seq2 = ''
    consensus_seq_by_anchor = []

    # print(seq_len)
    anchor_points = []
    anchor_number = 3
    for i in range(seq_len):
        try:
            seq_i = [ seq[i] for seq in list(sequences.values()) ]
        except IndexError:
            print(ids)
        # 第一类结果，不算gap，每个取出现概率较大的那个，并且概率要在80%以上
        # consensus_seq
        seq_i_A = seq_i.count('A')
        seq_i_T = seq_i.count('T')
        seq_i_C = seq_i.count('C')
        seq_i_G = seq_i.count('G')
        seq_i_GAP = seq_i.count('-')
        seq_i_nucle = sorted([('A', seq_i_A), ('T', seq_i_T), ('C', seq_i_C), ('G', seq_i_G), ('-', seq_i_GAP)], key = itemgetter(1), reverse=True)
        
        # test0
        # seq = seq_i_nucle[0][0]
        # seq2 = seq_i_nucle[0][0]

        # anchor point


        ### test2
        if i < 0  or i >= seq_len - 15:
            if seq_i_nucle[0][0] == '-':
               if seq_i_nucle[1][1] >= 2 and seq_i_nucle[1][1] > seq_i_nucle[2][1]:
                    consensus_seq_by_anchor.append(seq_i_nucle[1][0]) 
            else:
                consensus_seq_by_anchor.append(seq_i_nucle[0][0]) 
        
        if i >= 0 and i < seq_len - 15:
            #print(len(seq_i),seq_i_GAP,seq_i_nucle[0][1])
            if seq_i_nucle[0][0] == '-' :
                if seq_i_nucle[1][1] / seq_i_nucle[0][1] < 0.5:
                    consensus_seq_by_anchor.append(seq_i_nucle[0][0]) 
                else:
                    consensus_seq_by_anchor.append(seq_i_nucle[1][0]) 
                    
            else:
                if len(seq_i) - seq_i_GAP > seq_i_nucle[0][1]:
                    consensus_seq_by_anchor.append(seq_i_nucle[0][0]) 
                else:
                    anchor_point = i
                    anchor_points.append(anchor_point)
                    
                    if len(anchor_points) == 1:
                        consensus_seq_by_anchor.append(seq_i_nucle[0][0]) 
                        continue
                    elif len(anchor_points) <= anchor_number:
                        anchor_point = anchor_points[0]
                    else:
                        anchor_point = anchor_points[-anchor_number]
                    
                    seq_ij_max = get_seq_ij_count(seq_info, range(len(seq_info)), anchor_point, i )
                    # print(seq_ij_max)
                    
                    
                    # print('>>>')
                    gap_len =  ( i+1 - anchor_point) - len(seq_ij_max[0][0])
                    
                    if len(anchor_points) >=anchor_number:
                        pre_anchor_seq = "".join(consensus_seq_by_anchor[anchor_points[-anchor_number]:anchor_points[-2]+1]).replace('-','')
                    elif len(anchor_points) ==2:
                        pre_anchor_seq = "".join(consensus_seq_by_anchor[anchor_points[0]]).replace('-','')
                    #elif len(anchor_points) ==anchor_number-1:
                    #    pre_anchor_seq = "".join(consensus_seq_by_anchor[anchor_points[-3]:anchor_points[-2]+1]).replace('-','')

                    if seq_ij_max[0][0][0:len(pre_anchor_seq)] == pre_anchor_seq:
                        seq_expend = seq_ij_max[0][0][len(pre_anchor_seq):]
                    else:
                        # print("???")
                        # seq_ij_max1 = get_seq_ij_count(seq_info, range(len(seq_info)), anchor_point, i )
                        # print(seq_ij_max1)

                        anchor_point = anchor_points[-2]
                        seq_ij_max = get_seq_ij_count(seq_info, range(len(seq_info)), anchor_point, i )
                        seq_expend = seq_ij_max[0][0][1:]
                        
                    gap_len = anchor_points[-1] - anchor_points[-2] - len(seq_expend)
                    if gap_len > 0:
                        seq_expend = list('-'*gap_len + seq_expend)
                        # seq_ij = [seq_ij_temp[0],['-']*gap_len,seq_ij_temp[1:]]
                    else:
                        seq_expend = list(seq_expend)
                    if len(anchor_points) >=3:
                        consensus_seq_by_anchor = consensus_seq_by_anchor[:anchor_points[-2]+1]
                        consensus_seq_by_anchor.extend(seq_expend)
                    else:
                        consensus_seq_by_anchor = consensus_seq_by_anchor[:anchor_points[0]+1]
                        consensus_seq_by_anchor.extend(seq_expend)

                    # print(consensus_seq_by_anchor)
                    # print('\n')
                    pre_anchor_list = seq_ij_max[0][2]

                
        # print(consensus_seq_by_anchor)
        # print(anchor_points)


     
    # print("".join(consensus_seq_by_anchor))

    # print(consensus_seq)
    # print(ids)
    return [ids, "".join(consensus_seq_by_anchor).replace('-',''), "".join(consensus_seq)]
    # return [ids, consensus_seq2, "".join(consensus_seq)]
    # return [ids, "".join(consensus_seq), consensus_seq2]

def get_consensus_seq_mafft(align_seq_file, fre_cuoff, type, out_file, model):
    align_sequences = {}
    id = None
    seq = '' 
    aln_file = open('/home/boxu/temp/wtdbg2/for_qc/test.aln', 'w')
    for line in open(mafft( align_seq_file, model ), 'r') :
        if line.startswith('>'):
            if id != None:
                align_sequences[id] = seq.upper()
                aln_file.write(seq.upper()+"\n")
            id = line.strip().split('>')[1]
            seq = ''
        else:
            seq = seq + line.strip()
    align_sequences[id] = seq.upper()

    aln_file.write( seq.upper()+"\n")
    aln_file.close()

    if type == "seq":
        # align_seq = '\t|\t'.join(consensus_seq(align_sequences, fre_cuoff, type))
        align_seq_result = consensus_seq(align_sequences, fre_cuoff, type)
        align_seq = align_seq_result[1]
        g = open(out_file, 'w')
        g.write(">ctg1 "+str(len(align_seq)) + "\n" + align_seq)
        g.close()
        #aln_file = open('/home/boxu/temp/wtdbg2/for_qc/test.aln', 'a')
        #aln_file.write(align_seq_result[3] + "\n")
        #aln_file.close()
        #test_file.write('>'+align_seq_result[1]+'\n')
        #test_file.close()
    elif type == "tsd":
        # print(type)
        align_seq = consensus_seq(align_sequences, fre_cuoff,type)[1]
    return align_seq




def get_consensus_seq(consensus_method, consensus_seq_temp_file, out_path, ref_name, c_id, model ):
    # consensus_method = "mafft"
    if consensus_method == 'mafft':
        consensus_seq = get_consensus_seq_mafft( consensus_seq_temp_file, 0.75, 'seq',  out_path + ref_name + "/" + c_id + ".consensus.fa", model)
    elif consensus_method == 'wtdbg2':
        consensus_seq = get_consensus_seq_wtdbg2( consensus_seq_temp_file, out_path, ref_name, c_id, model )
    else:
        raise Exception("no this method")
    return consensus_seq



    