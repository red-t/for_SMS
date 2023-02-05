#!/usr/bin/python
import argparse
import pysam
import glob


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")
    parser.add_argument("--workdir", type=str, dest="workdir", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma", help="directory of the simulated dataset")
    parser.add_argument("--bam", type=str, dest="bam", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma/line_28_pacbio_gamma.bam", help="path to the alignment in BAM format.")
    parser.add_argument("--flank", type=int, dest="flank", default=50, help="flank size to extend around breakpoint upstream/downstream")
    parser.add_argument("--outdir", type=str, dest="outdir", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/test_with_simulation", help="output directory")
    return parser


def filter_high_freq_ins(sumry_f, min_freq=0.0001):
    ins_dict = {}
    for l in open(sumry_f, "r"):
        l = l.strip().split("\t")
        if float(l[6]) > min_freq:
            ins_dict[l[3]] = (l[0], int(l[1]), int(l[2]))

    return ins_dict


def find_target_hg(pgd_f, ins_dict, ins2hg, sub_size=50):
    order = int(pgd_f.split(".")[-2])
    for l in open(pgd_f, "r"):
        l = l.strip().split()
        if (l[0] != "#") and ("=" not in l[0]):
            for idx in range(1,sub_size+1):
                ins_id = l[idx]
                if (ins_id != "*") and (ins_id in ins_dict) :
                    hg = order*sub_size + int(idx)
                    if l[idx] in ins2hg:
                        ins2hg[ins_id].append(hg)
                    else:
                        ins2hg[ins_id] = [hg]


def filter_true_sup_aln(ins_dict, ins2hg, alnf, outdir="./", flank=50):
    for ins_id in ins_dict:
        ins = ins_dict[ins_id]
        true_outf_name = "{0}/true_{1}:{2}-{3}_{4}.bam".format(outdir, ins[0], ins[1], ins[2], ins_id)
        false_outf_name = "{0}/false_{1}:{2}-{3}_{4}.bam".format(outdir, ins[0], ins[1], ins[2], ins_id)
        true_outf = pysam.AlignmentFile(true_outf_name, "wb", template=alnf)
        false_outf = pysam.AlignmentFile(false_outf_name, "wb", template=alnf)
        for read in alnf.fetch(ins[0], ins[1], ins[2]):
            hg = read.query_name.split("hg")[1]
            hg = int(hg.split(":")[0])
            if hg in ins2hg[ins_id]:
                true_outf.write(read)
            else:
                false_outf.write(read)
        
        true_outf.close()
        false_outf.close()


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    alnf = pysam.AlignmentFile(args.bam, "rb")
    for chr in alnf.references:
        sumryf = "{0}/{1}/{2}.ins.summary".format(args.workdir, chr, chr)
        ins_dict = filter_high_freq_ins(sumryf, min_freq=0.0001)
        pgds = glob.glob("{0}/{1}/*pgd".format(args.workdir, chr))
        ins2hg = {}
        for pgd in pgds:
            find_target_hg(pgd, ins_dict, ins2hg, sub_size=50)
        
        filter_true_sup_aln(ins_dict, ins2hg, alnf, args.outdir, args.flank)

    alnf.close()


# python Filter_TrueSupportAlignment.py
