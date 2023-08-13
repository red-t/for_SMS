import argparse
from BPG_utils import build_popg

parser = argparse.ArgumentParser(description="")
parser.add_argument("--chassis", type=str, required=False, dest="ref_fasta", default=None, help="the chassis, i.e. the sequence into which TEs will be inserted; a fasta file")
parser.add_argument("--te-seqs", type=str, required=False, dest="te_fasta", default=None, help="TE sequences in a fasta file")
parser.add_argument("--pgd", type=str, required=True, dest="pgd_definition", default=None, help="the definition of the population genome")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file; will be multi-fasta file")
parser.add_argument("--ins-seq", type=str, required=False, dest="ins_seq", default=None, help="the output file of insertion sequence; will be multi-fasta file")
parser.add_argument("--sub_idx", type=int, required=True, dest="sub_idx", default=None, help="the index of the sub-population genome")
parser.add_argument("--sub_size", type=int, required=True, dest="sub_size", default=None, help="the size of the sub-population genome")
args = parser.parse_args()

build_popg(args.te_fasta,
           args.pgd_definition,
           args.ref_fasta,
           args.sub_idx,
           args.sub_size,
           args.output,
           args.ins_seq)
