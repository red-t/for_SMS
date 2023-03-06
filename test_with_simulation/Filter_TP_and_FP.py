#!/usr/bin/python
import argparse
import pysam
import glob
from subprocess import Popen, PIPE, STDOUT, DEVNULL


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")
    parser.add_argument("--workdir", type=str, dest="workdir", default="/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma", help="directory of the simulated dataset")
    parser.add_argument("--bam", type=str, dest="bam", default="./all_candidates_alignments.bam", help="path to the alignment in BAM format.")
    parser.add_argument("--qbed", type=str, dest="qbed", default="./all_candidates.bed", help="path to the insertion candidates BED file.")
    parser.add_argument("--rbed", type=str, dest="rbed", default="./AllIns.bed.gz", help="path to the simulated insertion information BED file, with tabix index.")
    parser.add_argument("--segf", type=str, dest="segf", default="./all_candidates_segments.txt", help="file of candidates segments information.")
    return parser

# ---------------------------------------------------------------

def classification(qbed, rbed):
    '''filter out tp/fp from all candidates
    
    Parameters
    ----------
    qbed: str
        all insertion candidates in bed format
    
    rbed: str
        all simulated insertion information in bgzip compressed bed
        format, with tabix index.

    Returns
    -------
    tp_dict: dict
        a dictionary for TP, insertion_ID <-- (chrom, start, end) 

    fp_dict: dict
        a dictionary for FP, insertion_ID <-- (chrom, start, end) 
    '''

    tabix        = pysam.TabixFile(rbed)
    tp_outf      = open("TP.bed", "w")
    fp_outf      = open("FP.bed", "w")
    special_outf = open("Special.bed", "w")

    tp_dict = {}
    fp_dict = {}
    for l in open(qbed, "r"):
        n = 0
        l = l.strip().split()

        if len(l):
            chrom = l[0]
            start = int(l[1]) - 50
            end   = int(l[2]) + 50
        else:
            break
        
        if start < 0:
            start = 0

        ite = tabix.fetch(chrom, start, end, parser=pysam.asTuple())
        while 1:
            try:
                row = ite.next()
                n  += 1
                tp_outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                              chrom, l[1], l[2], row[3], l[4], l[5], l[3]))
                tp_dict[l[3]] = set()
            except StopIteration:
                if n == 0:
                    fp_outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                chrom, l[1], l[2], l[3], l[4], '*'))
                    fp_dict[l[3]] = set()
                break
        
        if n>1:
            special_outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                               chrom, l[1], l[2], l[3], l[4], '*'))

    tp_outf.close()
    fp_outf.close()
    special_outf.close()
    return tp_dict, fp_dict

# ---------------------------------------------------------------

def tp_filter(segf, tp_dict, idxreads, rbf):

    # fetch all support alignments of tp candidates
    for l in open(segf, "r"):
        l = l.split()
        if l[7] in tp_dict:
            tp_dict[l[7]].add((l[16], int(l[11])))

    # filter & write out tp support alignments
    for id in tp_dict:

        # open output file for each tp candidate
        wbfn = 'tp_{}.bam'.format(id)
        wbf = pysam.AlignmentFile(wbfn, "wb", template=rbf)
        
        # filter tp support alignment based on it's segment information
        for qname, rpos in tp_dict[id]:
            ite = idxreads.find(qname)
            while 1:
                try:
                    aln = next(ite)
                    if aln.reference_start == rpos:
                        wbf.write(aln)
                except StopIteration:
                    break
        
        wbf.close()
        cmd_idx = ['samtools sort -o tmp.bam {} && mv tmp.bam {} && samtools index {}'.format(wbfn, wbfn, wbfn)]
        idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitcode = idx_proc.wait()

# ---------------------------------------------------------------

def fp_filter(segf, fp_dict, idxreads, rbf):

    # fetch all support alignments of fp candidates
    for l in open(segf, "r"):
        l = l.split()
        if l[7] in fp_dict:
            fp_dict[l[7]].add((l[16], int(l[11])))

    # filter & write out fp support alignments
    for id in fp_dict:

        # open output file for each fp candidate
        wbfn = 'fp_{}.bam'.format(id)
        wbf = pysam.AlignmentFile(wbfn, "wb", template=rbf)
        
        # filter fp support alignment based on it's segment information
        for qname, rpos in fp_dict[id]:
            ite = idxreads.find(qname)
            while 1:
                try:
                    aln = next(ite)
                    if aln.reference_start == rpos:
                        wbf.write(aln)
                except StopIteration:
                    break
        
        wbf.close()
        cmd_idx = ['samtools sort -o tmp.bam {} && mv tmp.bam {} && samtools index {}'.format(wbfn, wbfn, wbfn)]
        idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
        exitcode = idx_proc.wait()

# ---------------------------------------------------------------

def init_tp_dict(chr):
    '''initialize tp_dict from TP.bed

    Returns
    -------
    tp_dict: dict
        a dictionary for TP, insertion_ID <-- (chrom, start, end, 
        origin_id)
    '''
    tp_dict = {}
    for l in open("TP.bed", "r"):
        l = l.strip().split()
        if l[0] == chr:
            if l[3] in tp_dict:
                tp_dict[l[3]].append((l[0], int(l[1]), int(l[2]), l[6]))
            else:
                tp_dict[l[3]] = [(l[0], int(l[1]), int(l[2]), l[6])]

    return tp_dict

# ---------------------------------------------------------------

def find_target_hg(pgdf, tp_dict, ins2hg, sub_size=50):
    '''find target hg(s) for TP insertion
    
    Parameters
    ----------
    pgdf: str
        pgd (population genome definition) file, each pgd file cor-
        respond to sub_size hgs
    
    tp_dict: dict
        a dictionary for TP, insertion_ID <-- (chrom, start, end, 
        origin_id)
    
    ins2hg: dict
        a dictionary, insertion_ID <-- {hgx, hgy, hgz, ...}
    
    sub_size: int
        sub-population genome size, each pgdf contains information 
        of sub_size hgs

    Returns
    -------
    ins2hg: dict
        this function will be called many times, using different pgd-
        f to fill ins2hg. And return it as input of next calling.
    '''
    # pgdf comes from order-th sub-population genome
    order = int(pgdf.split(".")[-2])

    # read through body part of pgdf
    for l in open(pgdf, "r"):
        l = l.strip().split()
        if (l[0] != "#") and ("=" not in l[0]):
            for idx in range(1,sub_size+1):
                id = l[idx]
                if (id != "*") and (id in tp_dict):
                    hg = order*sub_size + int(idx)
                    if l[idx] in ins2hg:
                        ins2hg[id].add(hg)
                    else:
                        ins2hg[id] = set([hg])

# ---------------------------------------------------------------

def deep_filter(tp_dict, ins2hg):
    '''fiter tp alignments based on hg

    filter true support/unsupport alignments based on tp_dict a-
    nd ins2hg. The source hg of each alignment(read) can be infe-
    rred from alignment's qname. if alignment's source query name
    can match the target hgs of the insertion nearby, the alignm-
    ent would be considered as true support alignment.

    Parameters
    ----------
    tp_dict: dict
        tp_dict: dict
        a dictionary for TP, insertion_ID <-- (chrom, start, end, 
        origin_id)
    
    ins2hg: dict
        a dictionary, insertion_ID <-- {hgx, hgy, hgz, ...}
    '''
    for id in tp_dict:
        insl = tp_dict[id]
        for ins in insl:
            n   = 0
            oid = ins[3]

            # open file for reading
            fn  = 'tp_{}.bam'.format(oid)
            rbf = pysam.AlignmentFile(fn, 'rb', threads=5)
            
            # open file for writting
            tp_tfn  = 'tp_t_{}.bam'.format(oid)
            tp_ffn  = 'tp_f_{}.bam'.format(oid)
            tp_t    = pysam.AlignmentFile(tp_tfn, "wb", template=rbf)
            tp_f    = pysam.AlignmentFile(tp_ffn, "wb", template=rbf)

            # filter true/false support alignments
            # based on target & source hg
            for read in rbf.fetch():
                hg = read.query_name.split("hg")[1]
                hg = int(hg.split(":")[0])
                if hg in ins2hg[id]:
                    tp_t.write(read)
                else:
                    n += 1
                    tp_f.write(read)
            
            # close file
            tp_t.close()
            tp_f.close()
            rbf.close()

            # build index
            if n == 0:
                cmd_idx = ['mv tp_{}.bam tp_t_{}.bam && mv tp_{}.bam.bai tp_t_{}.bam.bai && rm tp_f_{}.bam'.format(oid, oid, oid, oid, oid)]
            else:
                cmd_idx = ['rm tp_{}.bam && rm tp_{}.bam.bai && samtools index tp_t_{}.bam && samtools index tp_f_{}.bam'.format(oid, oid, oid, oid)]

            idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
            exitcode = idx_proc.wait()
            if exitcode != 0:
                raise Exception("Error: samtools index for {} failed".format(oid))

# ---------------------------------------------------------------

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    
    # open candidates BAM file for reading
    rbf  = pysam.AlignmentFile(args.bam, "rb", threads=5)
    refs = rbf.references

    # build IndexedReads object
    idxreads = pysam.IndexedReads(rbf)
    idxreads.build()

    # classify TP/FP based on ground truth
    tp_dict, fp_dict = classification(args.qbed, args.rbed)

    # primary filter
    tp_filter(args.segf, tp_dict, idxreads, rbf)
    fp_filter(args.segf, fp_dict, idxreads, rbf)

    # close
    rbf.close(); del rbf
    
    # --------------------------------

    # process by chroms
    for chr in refs:
        # initialize tp_dict & ins2hg
        tp_dict = init_tp_dict(chr)
        ins2hg  = dict()
        if len(tp_dict) == 0:
            continue

        # fill ins2hg with pgd files
        pgdfs = glob.glob("{0}/{1}/*pgd".format(args.workdir, chr))
        for pgdf in pgdfs:
            find_target_hg(pgdf, tp_dict, ins2hg)
        
        # perform deep filteration
        deep_filter(tp_dict, ins2hg)


# python Filter_TP_and_FP.py