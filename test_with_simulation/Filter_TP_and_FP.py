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
    parser.add_argument("--subsize", type=int, dest="subsize", default=50, help="size of sub population genome")
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
    fp_outf      = open("FP.bed", "w")
    special_outf = open("Special.bed", "w")

    oid_id  = {}
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
        
        # for GRCh38 somatic insertions
        try:
            ite = tabix.fetch(chrom, start, end, parser=pysam.asTuple())
        except ValueError:
            fp_outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                          chrom, l[1], l[2], l[3], l[4], '*'))
            fp_dict[l[3]] = set()
            continue

        while 1:
            try:
                row = next(ite)
                n  += 1
                # oid <- [tchrom, tstart, tend, tid, nseg, tstrand, ostart, oend, oid, nconsist]
                oid_id[l[3]]  = [row[0], row[1], row[2], row[3], l[4], row[5], l[1], l[2], l[3], 0]
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

    fp_outf.close()
    special_outf.close()
    return tp_dict, fp_dict, oid_id

# ---------------------------------------------------------------

def tp_filter(segf, tp_dict, fp_dict, oid_id, idxreads, rbf):

    # fetch all support alignments of tp candidates
    for l in open(segf, "r"):
        l     = l.split()
        chr   = l[0]
        oid   = l[1]
        rpos  = int(l[2])
        qname = l[3]
        
        # determine if candidate of segment overlaps with ground truth
        if oid in tp_dict:
            tp_dict[oid].add((qname, rpos))
            # determine if chr of segment consist with candidate
            qchr  = qname.split("_")[0]
            qchr  = qchr.split(";")[1]
            if chr == qchr:
                oid_id[oid][-1] += 1
    # ---------------------------------------------------
    # open file for writting
    fp_outf = open("FP.bed", "a")
    tp_outf = open("TP.bed", "a")
    oids    = list(tp_dict.keys())
    for oid in oids:
        if oid_id[oid][-1] == 0:
            tp_dict.pop(oid)
            fp_outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                          oid_id[oid][0], oid_id[oid][6], oid_id[oid][7],
                          oid_id[oid][8], oid_id[oid][4], "*"))
            fp_dict[oid] = set()
        else:
            tp_outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                          oid_id[oid][0], oid_id[oid][1], oid_id[oid][2],
                          oid_id[oid][3], oid_id[oid][4], oid_id[oid][5],
                          oid_id[oid][6], oid_id[oid][7], oid_id[oid][8]))

    # close
    fp_outf.close()
    tp_outf.close()
    # ---------------------------------------------------
    # filter & write out tp support alignments
    for oid in tp_dict:

        # open output file for each tp candidate
        wbfn = 'tp_{}.bam'.format(oid)
        wbf  = pysam.AlignmentFile(wbfn, "wb", template=rbf)
        
        # filter tp support alignment based on it's segment information
        for qname, rpos in tp_dict[oid]:
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
        l     = l.split()
        oid   = l[1]
        rpos  = int(l[2])
        qname = l[3]
        if oid in fp_dict:
            fp_dict[oid].add((qname, rpos))

    # filter & write out fp support alignments
    for oid in fp_dict:

        # open output file for each fp candidate
        wbfn = 'fp_{}.bam'.format(oid)
        wbf = pysam.AlignmentFile(wbfn, "wb", template=rbf)
        
        # filter fp support alignment based on it's segment information
        for qname, rpos in fp_dict[oid]:
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

def init_tp_dict(chrom):
    '''initialize tp_dict from TP.bed

    Returns
    -------
    tp_dict: dict
        a dictionary for TP, insertion_ID <-- (chrom, start, end, 
        origin_id)
    '''
    tp_dict = {}
    for l in open("TP.bed", "r"):
        l     = l.strip().split()
        chr   = l[0]
        start = int(l[1])
        end   = int(l[2])
        id    = l[3]
        oid   = l[8]

        if chr == chrom:
            if id in tp_dict:
                tp_dict[id].append((chr, start, end, oid))
            else:
                tp_dict[id] = [(chr, start, end, oid)]

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

def deep_filter(tp_dict, ins2hg, chr):
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
    logf = open('log.txt', 'a')
    for id in tp_dict:
        insl = tp_dict[id]
        for ins in insl:
            n   = 0
            m   = 0
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
                # filter by qchr & chr
                qchr = read.query_name.split("_")[0]
                qchr = qchr.split(";")[1]
                if qchr != chr:
                    n += 1
                    tp_f.write(read)
                    continue
                # filter by source & target hg
                hg = read.query_name.split("hg")[1]
                hg = int(hg.split(":")[0])
                if hg in ins2hg[id]:
                    m += 1
                    tp_t.write(read)
                else:
                    n += 1
                    tp_f.write(read)
            
            # close file
            tp_t.close()
            tp_f.close()
            rbf.close()

            # build index
            if   n == 0:
                # print("all support alignments of {} is true".format(oid))
                cmd_idx = ['mv tp_{}.bam tp_t_{}.bam && mv tp_{}.bam.bai tp_t_{}.bam.bai && rm tp_f_{}.bam'.format(oid, oid, oid, oid, oid)]
            elif m == 0:
                print("all support alignments of {} is false".format(oid), file=logf)
                cmd_idx = ['mv tp_{}.bam fp_{}.bam && mv tp_{}.bam.bai fp_{}.bam.bai && rm tp_*_{}.bam'.format(oid, oid, oid, oid, oid)]
            else:
                # print("there are both true & flase support alignments of {}".format(oid))
                cmd_idx = ['rm tp_{}.bam && rm tp_{}.bam.bai && samtools index tp_t_{}.bam && samtools index tp_f_{}.bam'.format(oid, oid, oid, oid)]

            idx_proc = Popen(cmd_idx, stderr=DEVNULL, shell=True, executable='/bin/bash')
            exitcode = idx_proc.wait()
            if exitcode != 0:
                raise Exception("Error: samtools index for {} failed".format(oid))
    
    logf.close()

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
    tp_dict, fp_dict, oid_id = classification(args.qbed, args.rbed)

    # primary filter
    tp_filter(args.segf, tp_dict, fp_dict, oid_id, idxreads, rbf)
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
            find_target_hg(pgdf, tp_dict, ins2hg, args.subsize)
        
        # perform deep filteration
        deep_filter(tp_dict, ins2hg, chr)


# python Filter_TP_and_FP.py
