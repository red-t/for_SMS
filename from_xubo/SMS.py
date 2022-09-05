
from pysam import AlignmentFile
from concurrent.futures import ThreadPoolExecutor, as_completed
from read_alignment import build_cluster, process_cluster
from mappy import revcomp
import subprocess
import sys

bamFile=sys.argv[1]
out_path=sys.argv[2]
te_index=sys.argv[3]
te_size=sys.argv[4]
flanksize=sys.argv[5]
# te_size="/data/tusers/boxu/annotation/dm6_clean/dm6_clean.transposon.size"
# te_size="/data/tusers/boxu/annotation/hs37d5/hs37d5.transposon.size"
temp_out_path=out_path+"/temp/"

read_seq_dic = {}
te_size_dic = {}
for line in open(te_size, 'r'):
    line = line.strip().split('\t')
    te_size_dic[line[0]] = line[1]



def main(bam=bamFile, te_index=te_index):
    chrom2clusters = dict()
    with ThreadPoolExecutor(max_workers=5) as executor:
        bam_file = AlignmentFile(bam, 'rb')
        chroms = list(bam_file.references)
        # chroms = ['chr2L', 'chr3R']

        for read in AlignmentFile("/data/tusers/boxu/lrft/result/simulation/10/line_28_pacbio.10X.q0.sorted.bam", 'rb'):
        # for read in AlignmentFile(bamFile, 'rb'):
            if not read.is_supplementary and not read.is_secondary and read.query_sequence:
                # 以genome正向为标准
                # read_seq_dic[read.query_name] = read.query_sequence
                if read.is_reverse:
                    read_seq_dic[read.query_name] = revcomp(read.query_sequence)
                else:
                    read_seq_dic[read.query_name] = read.query_sequence
        # samtools根据read name提取alignment
        future2chrom = {executor.submit(build_cluster, bam, chrom, te_index, temp_out_path, read_seq_dic, te_size_dic, int(flanksize)):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            chrom2clusters[chrom] = future.result()
        
        future2chrom = {executor.submit(process_cluster, chrom2clusters[chrom], chrom, out_path):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            chrom2clusters[chrom] = future.result()
    merge_script = "cat " + out_path + "/*tmp.bed" + " > " + out_path + "/SMS.insertion.bed"
    subprocess.Popen(merge_script, shell=True)
    
    return chrom2clusters



chrom2c = main()
