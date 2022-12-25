

from pysam import AlignmentFile
from concurrent.futures import ThreadPoolExecutor, as_completed
from read_alignment import build_cluster, process_cluster
from mappy import revcomp
import subprocess
import sys
import sys
import argparse


def parse_args():
    """parse args for test"""

    parser = argparse.ArgumentParser(add_help=False)
    detect_insertion_setting = parser.add_argument_group('detect_insertion_setting')

    detect_insertion_setting.add_argument('-b', '--bam', dest='bamFile', type=str,
                                 help='Path of bam file, mapped by minimap2 -Y', default='')
    detect_insertion_setting.add_argument('-o', '--out_path', dest='out_path', type=str,
                                 help='Path of the output', default='./')
    detect_insertion_setting.add_argument('-t', '--te_index', dest='te_index', type=str, nargs='*',
                                 help='Transposon reference sequence index build by minimap2', default='')
    detect_insertion_setting.add_argument('-s', '--te_size', dest='te_size', type=str,
                                 help='Path of ransposon reference size file', default='')
    detect_insertion_setting.add_argument('-f', '--flanksize', dest='flanksize', type=str,
                                 help='Flanksize to the breakpoint', default='')
    detect_insertion_setting.add_argument('-p', '--prefix', dest='prefix', type=str,
                                 help='Path of ransposon reference size file', default='')
    detect_insertion_setting.add_argument('-g', '--genome_fa', dest='genome_fa', type=str,
                                 help='Path of ransposon reference size file', default='')                             
    detect_insertion_setting.add_argument('-i', '--genome_idx', dest='genome_idx', type=str,
                                 help='Path of ransposon reference size file', default='')
    detect_insertion_setting.add_argument('-a', '--te_anno_fa', dest='te_anno_fa', type=str,
                                 help='Path of ransposon reference size file', default='')
    detect_insertion_setting.add_argument('-r', '--repeatmasker_file', dest='repeatmasker_file', type=str,
                                 help='Path of ransposon reference size file', default='')
    detect_insertion_setting.add_argument('-h', '--help', dest='help', type=str,
                                 help='Help information', default='')


    return parser


def command_line_args(args):
    need_print_help = False if args else True
    parser = parse_args()
    args = parser.parse_args(args)

    if args.help or need_print_help:
        parser.print_help()
        sys.exit(1)
        
    return args





# check Q3: 有一些 module 没必要现在 import, 然后获取参数这部分其实可以用 argparse.ArgumentParser 来完善，具体可以看 simulation/my-define-landscape_random-insertions-freq-range.py


# check Q4: 下面编写的，名为 main 的函数，并不能起到防止错误调用的情形，之后如果要在不改变功能的前提下进行重构，应采用 if __name__ == '__main__' 这一判断
# 参考 zhihu.com/question/49136398
# simulation/my-define-landscape_random-insertions-freq-range.py 的结构可能比较相似
def main():

    bamFile=args.bamFile
    out_path=args.out_path
    te_index=args.te_index
    te_size=args.te_size
    flanksize=args.flanksize
    prefix=args.prefix
    genome_fa=args.genome_fa
    genome_idx=args.genome_idx
    te_anno_fa=args.te_anno_fa
    repeatmasker_file=args.repeatmasker_file

    temp_out_path=out_path+"/temp/"

    read_seq_dic = {}
    te_size_dic = {}
    for line in open(te_size, 'r'):
        line = line.strip().split('\t')
        te_size_dic[line[0]] = line[1]


    chrom2clusters = dict()
    with ThreadPoolExecutor(max_workers=5) as executor:
        bam_file = AlignmentFile(bamFile, 'rb')
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
        future2chrom = {executor.submit(build_cluster, bamFile ,genome_fa, chrom, te_index, temp_out_path, read_seq_dic, te_size_dic, int(flanksize)):chrom for chrom in chroms}
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


if __name__ == '__main__':
    args = command_line_args(sys.argv[1:])
    chrom2c = main()



