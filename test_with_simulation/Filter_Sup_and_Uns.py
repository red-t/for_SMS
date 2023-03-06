#!/usr/bin/python
import argparse
import pysam
import glob
from subprocess import Popen, PIPE, STDOUT, DEVNULL


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")
    parser.add_argument("--workdir", type=str, dest="workdir", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma", help="directory of the simulated dataset")
    parser.add_argument("--bam", type=str, dest="bam", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma/line_28_pacbio.bam", help="path to the alignment in BAM format.")
    parser.add_argument("--flank", type=int, dest="flank", default=50, help="flank size to extend around breakpoint upstream/downstream")
    parser.add_argument("--outdir", type=str, dest="outdir", default="./", help="output directory")
    return parser


def filter_high_freq_ins(sumry_f, min_freq=0.0001):
    '''filter insertion with frequency > min_freq from ground truth
    
    Parameters
    ----------
    sumry_f: str
        summary file of the simulated insertions
    
    min_freq: float
        min frequency, insertion with frequency > min_freq would be
        considered as "high frequency", in other words, "germline"

    Returns
    -------
    ins_dict: dict
        a dictionary, insertion_ID <-- (chrom, start, end)
    '''

    ins_dict = {}
    for l in open(sumry_f, "r"):
        l = l.strip().split("\t")
        if float(l[6]) > min_freq:
            ins_dict[l[3]] = (l[0], int(l[1]), int(l[2]))

    return ins_dict


def find_target_hg(pgd_f, ins_dict, ins2hg, sub_size=50):
    '''find target hg(s) for each insertion
    
    Parameters
    ----------
    pgd_f: str
        pgd (population genome definition) file, each pgd file cor-
        respond to sub_size hgs
    
    ins_dict: dict
        a dictionary, insertion_ID <-- (chrom, start, end)
    
    ins2hg: dict
        a dictionary, insertion_ID <-- {hgx, hgy, hgz, ...}
    
    sub_size: int
        sub-population genome size, each pgd_f contains information 
        of sub_size hgs

    Returns
    -------
    ins2hg: dict
        this function will be called many times, using different pgd-
        f to fill ins2hg. And return it as input of next calling.
    '''
    order = int(pgd_f.split(".")[-2])
    for l in open(pgd_f, "r"):
        l = l.strip().split()
        if (l[0] != "#") and ("=" not in l[0]):
            for idx in range(1,sub_size+1):
                ins_id = l[idx]
                if (ins_id != "*") and (ins_id in ins_dict) :
                    hg = order*sub_size + int(idx)
                    if l[idx] in ins2hg:
                        ins2hg[ins_id].add(hg)
                    else:
                        ins2hg[ins_id] = set([hg])


def filter_true_sup_aln(ins_dict, ins2hg, alnf, outdir="./", flank=50):
    '''filter true support alignment

    filter true support/unsupport alignments based on ins_dict a-
    nd ins2hg. The source hg of each alignment(read) can be infe-
    rred from alignment's qname. if alignment's source query name
    can match the target hgs of the insertion nearby, the alignm-
    ent would be considered as true support alignment.

    Parameters
    ----------
    ins_dict: dict
        a dictionary, insertion_ID <-- (chrom, start, end)
    
    ins2hg: dict
        a dictionary, insertion_ID <-- {hgx, hgy, hgz, ...}
    
    alnf: AlignmentFile
        pysam.AlignmentFile instance
    
    outdir: str
        output directory
    
    flank: int
        flank size, extend flank(bp) upstream and downstream the 
        breakpoints when fetching alignments
    '''
    for ins_id in ins_dict:
        ins     = ins_dict[ins_id]
        chrom   = ins[0]
        start   = ins[1] - flank
        end     = ins[2] + flank
        n       = 0
        m       = 0

        # open file for writting
        t_fn    = "{}/sup_{}_{}.bam".format(outdir, chrom, ins_id)
        f_fn    = "{}/uns_{}_{}.bam".format(outdir, chrom, ins_id)
        tf      = pysam.AlignmentFile(t_fn, "wb", template=alnf)
        ff      = pysam.AlignmentFile(f_fn, "wb", template=alnf)

        # filter sup/unsup alignments based on source & target hg
        for read in alnf.fetch(chrom, start, end):
            hg = read.query_name.split("hg")[1]
            hg = int(hg.split(":")[0])
            if hg in ins2hg[ins_id]:
                tf.write(read)
                n += 1
            else:
                ff.write(read)
                m += 1
        
        # close
        tf.close()
        ff.close()
        
        # index
        if n > 0:
            cmd_idx = ['samtools index {}'.format(t_fn)]
            idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
            exitcode = idx_proc.wait()
        if m > 0:
            cmd_idx = ['samtools index {}'.format(f_fn)]
            idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
            exitcode = idx_proc.wait()


if __name__ == '__main__':
    parser = get_parser()
    args   = parser.parse_args()
    alnf   = pysam.AlignmentFile(args.bam, "rb")
    refs   = alnf.references

    for chr in refs:
        sumryf   = "{0}/{1}/{2}.ins.summary".format(args.workdir, chr, chr)
        ins_dict = filter_high_freq_ins(sumryf, min_freq=0.0001)
        pgds     = glob.glob("{0}/{1}/*pgd".format(args.workdir, chr))
        ins2hg   = {}
        for pgd in pgds:
            find_target_hg(pgd, ins_dict, ins2hg, sub_size=50)
        
        filter_true_sup_aln(ins_dict, ins2hg, alnf, args.outdir, args.flank)

    alnf.close()


# python Filter_TrueSupportAlignment.py
