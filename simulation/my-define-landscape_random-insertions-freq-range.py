#!/usr/bin/python
import argparse
import random
import re


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument("--chassis", type=str, required=True, dest="ref_fa", default=None, help="the chassis, i.e. the sequence into which TEs will be inserted; a fasta file")
    parser.add_argument("--te-seqs", type=str, required=True, dest="te_fa", default=None, help="TE CCS sequences in a fasta file")
    parser.add_argument("--N", type=int, required=True, dest="pop_size", default=None, help="the number of haploid genomes")
    parser.add_argument("--divergence-rate", required=True, dest="d_rate", help="TE sequence divergence rate (percent)")
    parser.add_argument("--germline-count", type=int, required=True, dest="germline_count", default=1, help="germline insertions mumber in the population genome")
    parser.add_argument("--somatic-count", type=int, required=True, dest="somatic_count", default=1, help="somatic insertions mumber in the population genome")
    parser.add_argument("--output", type=str, required=False, dest="output", default=None, help="the output file")
    parser.add_argument("--min-distance", type=int, required=False, dest="min_distance", default=1000, help="minimum distance between TE insertions")
    
    return parser


def define_nest_ins(idx2te, parent_te_len, parent_trunc_len, te_type, d_rate):
    if te_type == "s":
        ccs_idx = random.choice([72, 15, 54, 44, 100, 106])
    else:
        ccs_idx = random.randint(1,len(idx2te))
    te_id = idx2te[ccs_idx][0]
    te_len = idx2te[ccs_idx][1]

    # TSD
    tsd = random.randint(6,8)

    # Position
    pos = random.randint(tsd, parent_te_len-parent_trunc_len)

    # Strand
    if random.random() < 0.5:
        strand = "-"
        l_strand = "R"
    else:
        strand = "+"
        l_strand = "F"

    # Truncations
    trunc_id = "U"
    trunc_len = 0
    trunc_dsl = ""
    if random.random() <= 0.01:
        trunc_id = "T"
        trunc_len = int(random.uniform(0.1, 0.6)*te_len)
        trunc_start = random.randint(1,te_len-trunc_len)
        trunc_dsl = "[{0}..{1}]".format(trunc_start, trunc_start + trunc_len - 1)

    # Divergence
    if d_rate > 0:
        # div_rate = round(random.uniform(0,d_rate), 2)
        div_rate = d_rate

    # Nest
    nest_dsl = ""
    nest_id = ""
    if random.random() <= 0.01:
        nest_id, nest_dsl=define_nest_ins(idx2te, te_len, trunc_len, te_type, d_rate)
        nest_id = "(" + nest_id + ")"
        nest_dsl = "{" + nest_dsl + "}"

    left = "{0}~{1}{2}{3}".format(te_id, trunc_id, l_strand, nest_id)
    if d_rate > 0:
        right = "{0}:${1}{2}{3}{4}{5}%{6}bp".format(pos, ccs_idx, trunc_dsl, strand, nest_dsl, div_rate, tsd)
    else:
        right = "{0}:${1}{2}{3}{4}{5}bp".format(pos, ccs_idx, trunc_dsl, strand, nest_dsl, tsd)
    return left, right


def define_ins(idx, idx2te, te_type="g", d_rate=0):
    if te_type == "s":
        ccs_idx = random.choice([72, 15, 54, 44, 100, 106])
    else:
        ccs_idx = random.randint(1,len(idx2te))
    te_id = idx2te[ccs_idx][0]
    te_len = idx2te[ccs_idx][1]

    # TSD
    tsd = random.randint(6,8)

    # Strand
    if random.random() < 0.5:
        strand = "-"
        l_strand = "R"
    else:
        strand = "+"
        l_strand = "F"

    # Truncations
    trunc_id = "U"
    trunc_len = 0
    trunc_dsl = ""
    if random.random() <= 0.1:
        trunc_id = "T"
        trunc_len = int(random.uniform(0.1, 0.6)*te_len)
        trunc_start = random.randint(1,te_len-trunc_len)
        trunc_dsl = "[{0}..{1}]".format(trunc_start, trunc_start + trunc_len - 1)

    # Divergence
    if d_rate > 0:
        # div_rate = round(random.uniform(0,d_rate), 2)
        div_rate = d_rate

    # Nest
    nest_dsl = ""
    nest_id = ""
    if random.random() <= 0.1:
        nest_id, nest_dsl=define_nest_ins(idx2te, te_len, trunc_len, te_type, d_rate)
        nest_id = "(" + nest_id + ")"
        nest_dsl = "{" + nest_dsl + "}"

    left = "{0}~{1}~{2}{3}{4}".format(idx, te_id, trunc_id, l_strand, nest_id)
    if d_rate > 0:
        right = "${0}{1}{2}{3}{4}%{5}bp".format(ccs_idx, trunc_dsl, strand, nest_dsl, div_rate, tsd)
    else:
        right = "${0}{1}{2}{3}{4}bp".format(ccs_idx, trunc_dsl, strand, nest_dsl, tsd)
    # te_dsl = "{0}={1}".format(left, right)
    return left, right


def define_header(ref_fa, te_fa, germline_count, somatic_count, germline_pos, d_rate):
    ref_idx=ref_fa+".fai"
    te_idx=te_fa+".fai"
    idx2te={}
    id2dsl={}

    # reference information
    for l in open(ref_idx, "r"):
        l=l.split()
        contig_id = l[0]
        genome_size = int(l[1])
        with open("{0}.tmp.pgd.header".format(contig_id), "w") as fout:
            fout.write("# Chasis {0}; Length {1} nt\n".format(contig_id, l[1]))
    
    # parse TE information
    counter=1
    for l in open(te_idx, "r"):
        l=l.split()
        if l[0] != "":
            idx2te[counter]=(l[0], int(l[1]))
            counter += 1

    # simualte total (germline_count + somatic_count) insertions
    ins_ids = []
    ins_count = germline_count + somatic_count
    fout = open("{0}.tmp.pgd.header".format(contig_id), "a")
    for idx in range(1, ins_count+1):
        if (idx-1) in germline_pos:
            left, right = define_ins(idx, idx2te, d_rate=d_rate)
        else:
            left, right = define_ins(idx, idx2te, "s", d_rate=d_rate)
        
        ins_ids.append(left)
        id2dsl[left] = right
        fout.write("{0}={1}\n".format(left, right))
    
    fout.close()
    return id2dsl, contig_id, genome_size, ins_ids


def define_body(id2dsl, popsize, germline_count, somatic_count, contig_id, mindist, maxdist, germline_pos, ins_ids):
    fout = open("{0}.tmp.pgd.body".format(contig_id), "w")
    summary = open("{0}.ins.summary".format(contig_id), "w")
    counter = 1

    # write insertions
    for i in range(germline_count+somatic_count):
        dist = random.randint(mindist,maxdist)
        counter += dist # in other words, "pos"
        if i in germline_pos:
            popfreq=random.uniform(0.1, 1)
            count_te = int(popfreq*popsize)
        else:
            popfreq = float(1)/popsize
            count_te = 1
        count_empty = popsize - count_te
        id = ins_ids[i]
        toshuf = [id for i in range(count_te)]
        toshuf.extend("*"*count_empty)
        random.shuffle(toshuf)
        toshuf.insert(0, str(counter))
        tows = " ".join(toshuf)
        # write out body
        fout.write(tows+"\n")

        # write out insertion information
        dsl = id2dsl[id]
        strand = re.search(r"([+-])", dsl).group(1)
        tsd_len = int(re.search("([\d\.]+)bp$", dsl).group(1))
        tsd_start = counter - tsd_len + 1
        summary.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(contig_id, counter, counter+1, id, dsl, strand, popfreq, tsd_start, counter))
    
    fout.close()
    summary.close()





if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    # print("my-define-landscape_random-insertions-freq-range.py parameters:")
    # print("ref_fa", args.ref_fa)
    # print("te_fa", args.te_fa)
    # print("pop_size", args.pop_size)
    # print("germline_count", args.germline_count)
    # print("somatic_count", args.somatic_count)
    # print("output", args.output)
    # print("min-distance", args.min_distance)
    if (args.germline_count + args.somatic_count) > 0:

        germline_pos = []
        while len(germline_pos) < args.germline_count:
            x = random.randint(0, args.germline_count + args.somatic_count)
            if x not in germline_pos:
                germline_pos.append(x)

        # define header of the pgd-file (population genome definition)
        id2dsl, contig_id, genomesize, ins_ids = define_header(args.ref_fa, args.te_fa, args.germline_count, args.somatic_count, germline_pos, float(args.d_rate))
        
        # define body of the pgd-file (population genome definition)
        mindist=args.min_distance
        maxdist=int(genomesize/(args.germline_count+args.somatic_count+1))
        if(maxdist<=mindist):
            raise ValueError("Genome of "+str(genomesize)+ " is too small for "+str(args.germline_count + args.somatic_count)+ " insertions with a min-distance between insertions of "+str(mindist))
        
        define_body(id2dsl, args.pop_size, args.germline_count, args.somatic_count, contig_id, mindist, maxdist, germline_pos, ins_ids)
