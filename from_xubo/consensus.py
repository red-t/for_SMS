import os
import subprocess
from operator import itemgetter
from collections import Counter
from subprocess import Popen, PIPE, STDOUT, DEVNULL
from pysam import  FastaFile



###          ###
###  wtdbg2  ###
###          ###
def get_consensus_seq_wtdbg2(consensus_seq_temp_file,out_path,ref_name, c_id):
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
    
    wtdbg2_consensus_cmd = [ 'wtpoa-cns', '-q', '-t', '10', '-i', out_path + ref_name + "/"+c_id +".wtdbg2.consensus.temp2.ctg.lay.gz", '-fo', out_path + ref_name + "/"+c_id +".consensus.fa" ]
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


###          ###
###   canu   ###
###          ###
def get_consensus_seq_canu(consensus_seq_temp_file, out_path, ref_name):
    canu_contig_cmd = [  'canu -p canu -fast correctedErrorRate=0.205 minInputCoverage=1 stopOnLowCoverage=1 minReadLength=200 minOverlapLength=100 corOutCoverage=1 -d', out_path + ref_name, 'genomeSize=8k', '-pacbio', consensus_seq_temp_file ]        
    canu_contig_proc = Popen([" ".join(canu_contig_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = canu_contig_proc.wait()

    rename_file_cmd = ['mv', out_path + ref_name + "/canu.contigs.fasta", out_path + ref_name+ "/insert.consensus.fa"]
    rename_file_proc = Popen([" ".join(rename_file_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash')
    exitcode = rename_file_proc.wait()

    # print(len(open(out_path + ref_name + ".consensus.temp.canu_consensus.fa", 'r').readlines()))
    if len(open(out_path + ref_name + "/canu.contigs.fasta", 'r').readlines()) == 0:
        consensus_seq = ''
    else:
        # fa_index
        # fa_index_cmd = ['samtools', 'faidx', out_path + ref_name + ".consensus.temp.canu_consensus.fa" ]
        # fa_index_proc = Popen([" ".join(fa_index_cmd)], stderr=DEVNULL, shell=True, executable='/bin/bash' )
        consensus_seq_fa = FastaFile(out_path + ref_name + "/canu.contigs.fasta")
        consensus_seq = consensus_seq_fa.fetch('ctg1')
    return consensus_seq



###           ###
###   mafft   ###
###           ###
def mafft(align_seq_file):
    ''' use MAFFT to create MSA '''

    out_file_msa = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
    
    # 根据序列的长度进行参数的选择
    #     if
    # args = [ 'mafft', '--randomseed', '1', align_seq_file ]
    # args = [ 'mafft', '--retree', '1', align_seq_file ]
    # args = [ 'mafft', '--retree', '2', align_seq_file ]
    #args = [ 'mafft', '--localpair', '--maxiterate', '1000', align_seq_file ]
    args = [ 'mafft', '--localpair', '--maxiterate', '1000', '--op', '1', align_seq_file ]
    # args = [ 'mafft', '--genafpair', '--maxiterate','1000', align_seq_file ]

    # args = [ 'mafft', '--retree', '2', '--op', '2',align_seq_file ]
 
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = FNULL)

    with open(out_file_msa, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            # if line[0] != '>':
            out_fa.write(line)
    return out_file_msa
   
def consensus_seq(sequences, fre_cuoff, type):
    ids = '*'
    if list(sequences.keys()) != [None]:
        ids = ';'.join(list(sequences.keys()))
    seq_len = len(list(sequences.values())[0])
    consensus_seq = []
    consensus_seq2 = ''

    # print(seq_len)

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
        seq = seq_i_nucle[0][0]
        seq2 = seq_i_nucle[0][0]

        # if seq_i_nucle[0][0] == '-':
         #   if seq_i_nucle[1][1] >= 2 and seq_i_nucle[1][1] > seq_i_nucle[2][1]:
        #        seq2 = seq_i_nucle[1][0]
            


        ### test2
        if i <= 15 or i >= seq_len - 15:
            if seq_i_nucle[0][0] == '-':
               if seq_i_nucle[1][1] >= 2 and seq_i_nucle[1][1] > seq_i_nucle[2][1]:
                    seq = seq_i_nucle[1][0]
                    seq2 = seq_i_nucle[1][0]

        if i > 15 and i < seq_len - 15:
            if seq_i_nucle[0][0] == '-' or ( seq_i_nucle[0][1] == seq_i_nucle[1][1] and seq_i_nucle[0][1] >= 2 )or ( seq_i_nucle[0][1] == seq_i_nucle[1][1] and seq_i_nucle[0][1] == seq_i_nucle[2][1] and seq_i_nucle[0][1] >= 2) : #) or seq_i_GAP >= 4 :
                # seq = seq_i_nucle[0][0] + seq_i_nucle[1][0]
                seq_temp_list = []
                for seq in list(sequences.values()):
                    l = 1;k = 1
                    seq_temp = seq[i]
                    while l < 10:
                        l = l + 1 
                        if seq[i-l] != '-':
                            if k <= 2 :
                                seq_temp = seq[i-l] + seq_temp
                                k = k + 1
                            else:
                                break
                    if len(seq_temp) < 3:
                        seq_temp = ''.join(['-']*(3-len(seq_temp))) + seq_temp

                    l = 1;k = 1
                    while l < 10:
                        l = l + 1 
                        if seq[i+l] != '-':
                            if k <= 2 :
                                seq_temp =  seq_temp + seq[i+l]
                                k = k + 1
                            else:
                                break
                    if len(seq_temp) < 5:
                        seq_temp =  seq_temp + ''.join(['-']*(5-len(seq_temp)))
                    seq_temp_list.append(seq_temp)

                # print(seq_temp_list)
                seq_ij_dic = Counter(seq_temp_list)
                seq_ij_max = sorted([(k,seq_ij_dic[k]) for k in seq_ij_dic], key = itemgetter(1), reverse=True)
                if seq_ij_max[0][0] == '-----':
                    seq = seq_ij_max[1][0] + "_" + seq_i_nucle[0][0] + seq_i_nucle[1][0]
                    seq2 = seq_ij_max[1][0][2]
                else:
                    seq = seq_ij_max[0][0] + "_" + seq_i_nucle[0][0] + seq_i_nucle[1][0]
                    seq2 = seq_ij_max[0][0][2]
                # print(i)
                # print(seq)
                # seq2 = seq_ij_max[0][0][2]
                #if seq == '-':
                #    seq = ''

        
        if seq2 == '-':
            seq2 = ''
            consensus_seq2 = consensus_seq2 + seq2
        else:
            consensus_seq2 = consensus_seq2 + seq2
        consensus_seq.append(seq)


        

    # print(consensus_seq)
    # print(ids)
    return [ids, consensus_seq2, "".join(consensus_seq)]
    # return [ids, "".join(consensus_seq), consensus_seq2]

def get_consensus_seq_mafft(align_seq_file, fre_cuoff, type, out_file):
    align_sequences = {}
    id = None
    seq = '' 
    test_file = open('/home/boxu/temp/wtdbg2/for_qc/test.txt', 'w')
    for line in open(mafft( align_seq_file ), 'r') :
        if line.startswith('>'):
            if id != None:
                align_sequences[id] = seq.upper()
                test_file.write(seq.upper()+"\n")
            id = line.strip().split('>')[1]
            seq = ''
        else:
            seq = seq + line.strip()
    align_sequences[id] = seq.upper()

    test_file.write(seq.upper()+"\n")
    test_file.close()

    if type == "seq":
        # align_seq = '\t|\t'.join(consensus_seq(align_sequences, fre_cuoff, type))
        align_seq_result = consensus_seq(align_sequences, fre_cuoff, type)
        align_seq = align_seq_result[1]
        g = open(out_file, 'w')
        g.write(">ctg1 "+str(len(align_seq)) + "\n" + align_seq)
        g.close()
        #test_file = open('/home/boxu/temp/wtdbg2/for_qc/test.txt', 'a')
        #test_file.write('>'+align_seq_result[1]+'\n')
        #test_file.close()
    elif type == "tsd":
        # print(type)
        align_seq = consensus_seq(align_sequences, fre_cuoff,type)[1]
    return align_seq




def get_consensus_seq(consensus_method, consensus_seq_temp_file, out_path, ref_name, c_id ):
    if consensus_method == 'mafft':
        consensus_seq = get_consensus_seq_mafft( consensus_seq_temp_file, 0.75, 'seq',  out_path + ref_name + "/" + c_id + ".consensus.fa")
    elif consensus_method == 'canu':
        consensus_seq = get_consensus_seq_canu( consensus_seq_temp_file, out_path, ref_name )
    elif consensus_method == 'wtdbg2':
        consensus_seq = get_consensus_seq_wtdbg2( consensus_seq_temp_file, out_path, ref_name, c_id )
    else:
        raise Exception("no this method")
    return consensus_seq



    
