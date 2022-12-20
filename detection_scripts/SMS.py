

from pysam import AlignmentFile
from concurrent.futures import ThreadPoolExecutor, as_completed
from read_alignment import build_cluster, process_cluster
from mappy import revcomp
import subprocess
import sys

# Q3: 有一些 module 没必要现在 import, 然后获取参数这部分其实可以用 argparse.ArgumentParser 来完善，具体可以看 simulation/my-define-landscape_random-insertions-freq-range.py
bamFile=sys.argv[1]
out_path=sys.argv[2]
te_index=sys.argv[3]
te_size=sys.argv[4]
flanksize=sys.argv[5]
prefix=sys.argv[6]
genome_fa=sys.argv[7]
genome_idx=sys.argv[8]
te_anno_fa=sys.argv[9]
repeatmasker_file=sys.argv[10]

temp_out_path=out_path+"/temp/"

read_seq_dic = {}
te_size_dic = {}
for line in open(te_size, 'r'):
    line = line.strip().split('\t')
    te_size_dic[line[0]] = line[1]


# Q4: 下面编写的，名为 main 的函数，并不能起到防止错误调用的情形，之后如果要在不改变功能的前提下进行重构，应采用 if __name__ == '__main__' 这一判断
# 参考 zhihu.com/question/49136398
# simulation/my-define-landscape_random-insertions-freq-range.py 的结构可能比较相似
def main(bam=bamFile, te_index=te_index):
    chrom2clusters = dict()
    with ThreadPoolExecutor(max_workers=5) as executor:
        bam_file = AlignmentFile(bam, 'rb')
        chroms = list(bam_file.references)
        # chroms = ['chr2L', 'chr3R']

        # 如果minimap2 的比对模式中没有 -Y，那就需要便利一遍先提取hard clip的信息，如果有的话就不需要了
        # for read in AlignmentFile("/data/tusers/boxu/lrft2/result/50/line_28.genome.50X.q0.sorted.bam", 'rb'):
        # # for read in AlignmentFile(bamFile, 'rb'):
        #     if not read.is_supplementary and not read.is_secondary and read.query_sequence:
        #         # 以genome正向为标准
        #         # read_seq_dic[read.query_name] = read.query_sequence
        #         if read.is_reverse:
        #             read_seq_dic[read.query_name] = revcomp(read.query_sequence)
        #         else:
        #             read_seq_dic[read.query_name] = read.query_sequence
        
        
        
        # samtools根据read name提取alignment
        future2chrom = {executor.submit(build_cluster, bam,genome_fa, chrom, te_index, temp_out_path, read_seq_dic, te_size_dic, int(flanksize)):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            chrom2clusters[chrom] = future.result()
        
        # Q5: 之后我们可以改成以 cluster 为单位进行并行，而不是以 chromosome 为单位？
        future2chrom = {executor.submit(process_cluster, chrom2clusters[chrom], chrom, out_path, genome_fa, genome_idx, te_anno_fa, repeatmasker_file, te_size_dic ):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            chrom2clusters[chrom] = future.result()
    merge_script = "cat " + out_path + "/*tmp.bed" + " > " + out_path + "/" + prefix + ".insertion.bed"
    subprocess.Popen(merge_script, shell=True)
    
    return chrom2clusters



chrom2c = main()


# correspond to Q10 & Q11 in read_alignment.pyx
a=[]
for i in range(10):
    if i >= 2:
        if i >= 5:
            continue
        pass

    if i >= 2:
        a.append(i)

# a == [2, 3, 4]