import argparse
from define_pgdf_utils import get_germline_pos, define_header, define_body


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
        germline_pos = get_germline_pos(args.germline_count, args.somatic_count)

        # define header of the pgd-file (population genome definition)
        id2dsl, contig_id, genomesize, ins_ids = define_header(args.ref_fa, args.te_fa, args.germline_count, args.somatic_count, germline_pos, float(args.d_rate))
        
        # define body of the pgd-file (population genome definition)
        mindist = args.min_distance
        maxdist = int(genomesize / (args.germline_count + args.somatic_count + 1))
        if(maxdist <= mindist):
            raise ValueError("Genome of "+str(genomesize)+ " is too small for "+str(args.germline_count + args.somatic_count)+ " insertions with a min-distance between insertions of "+str(mindist))
        
        ntotal = args.germline_count + args.somatic_count
        define_body(id2dsl, args.pop_size, ntotal, contig_id, mindist, maxdist, germline_pos, ins_ids)