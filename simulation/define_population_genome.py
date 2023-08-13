import argparse
from define_pgdf_utils import get_germline_pos, define_header, define_body


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument("--chassis", type=str, dest="ref_fa", default=None, help="the chassis, i.e. the sequence into which TEs will be inserted; a fasta file")
    parser.add_argument("--te-seqs", type=str, dest="te_fa", default=None, help="TE CCS sequences in a fasta file")
    parser.add_argument("--N", type=int, dest="pop_size", default=None, help="the number of haploid genomes")
    parser.add_argument("--divergence-rate", dest="d_rate", help="TE sequence divergence rate (percent)")
    parser.add_argument("--germline-count", type=int, dest="germline_count", default=1, help="germline insertions mumber in the population genome")
    parser.add_argument("--somatic-count", type=int, dest="somatic_count", default=1, help="somatic insertions mumber in the population genome")
    parser.add_argument("--output", type=str, dest="output", default=None, help="the output file")
    parser.add_argument("--min-distance", type=int, dest="min_distance", default=1000, help="minimum distance between TE insertions")
    parser.add_argument("--species", type=str, dest="species", default="human", help="the transposon library to the corresponding species will be used")
    
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if (args.germline_count + args.somatic_count) > 0:        
        germline_pos = get_germline_pos(args.germline_count, args.somatic_count)

        # define header of the pgd-file (population genome definition)
        id2dsl, contig_id, genomesize, ins_ids = define_header(args.ref_fa, args.te_fa, args.germline_count, args.somatic_count, germline_pos, float(args.d_rate), args.species)
        
        # define body of the pgd-file (population genome definition)
        mindist = args.min_distance
        try:
            maxdist = int(genomesize / (args.germline_count + args.somatic_count))
        except:
            maxdist = int(genomesize / (args.germline_count + args.somatic_count + 1))
        
        if(maxdist <= mindist):
            raise ValueError("Genome of "+str(genomesize)+ " is too small for "+str(args.germline_count + args.somatic_count)+ " insertions with a min-distance between insertions of "+str(mindist))

        ntotal = args.germline_count + args.somatic_count
        define_body(id2dsl, args.pop_size, ntotal, contig_id, mindist, maxdist, germline_pos, ins_ids, args.species)