from TGS_utils import generate_TGS
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--pg", type=str, required=True, dest="pop_gen", default=None, help="the population genome - a fasta file")
parser.add_argument("--read-length", type=int, required=False, dest="read_length", default=None, help="the mean read length")
parser.add_argument("--tgs-alpha", type=float, required=False, dest="tgs_alpha", default=None, help="alpha of gamma model applied for TGS length generation.")
parser.add_argument("--tgs-loc", type=float, required=False, dest="tgs_loc", default=None, help="loc of gamma model applied for TGS length generation.")
parser.add_argument("--tgs-beta", type=float, required=False, dest="tgs_beta", default=None, help="beta of gamma model applied for TGS length generation.")
parser.add_argument("--rld-file", type=str, required=False, dest="rldfile", default="", help="read length distribution file; will override --read-length and --std-dev")
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="the error rate of the reads; indels")
parser.add_argument("--error-fraction", type=str, required=False, dest="err_frac", help="the fraction of each type of errors, mis:del:ins")
parser.add_argument("--reads", type=int, required=True, dest="reads", default=None, help="the total number of reads")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="output - a fasta file")
parser.add_argument("--tgs-maxl", type=int, required=True, dest="tgs_maxl", default=None, help="Max length of TGS reads")
parser.add_argument("--tgs-minl", type=int, required=True, dest="tgs_minl", default=None, help="Min length of TGS reads.")
args = parser.parse_args()


print("Simulating PacBio Pool-Seq reads")
print("Reading the length of the population genome")
print("Getting mutator for PacBio reads with error rate {0} and error fraction {1}".format(args.error_rate, args.err_frac))

generate_TGS(args.pop_gen,
             args.error_rate,
             args.err_frac,
             args.reads,
             args.read_length,
             args.tgs_alpha,
             args.tgs_loc,
             args.tgs_beta,
             args.rldfile,
             args.fasta,
             args.tgs_minl,
             args.tgs_maxl)

print("Finished\n")

