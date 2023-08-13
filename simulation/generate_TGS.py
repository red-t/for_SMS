from TGS_utils import generate_PACBIO, generate_ONT
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--pg", type=str, required=True, dest="pop_gen", default=None, help="the population genome - a fasta file")
parser.add_argument("--reads", type=int, required=True, dest="nread", default=None, help="the total number of reads")
parser.add_argument("--fasta", type=str, required=True, dest="outfa", default=None, help="output - a fasta file")
parser.add_argument("--tgs-maxl", type=int, required=True, dest="maxl", default=None, help="Max length of TGS reads")
parser.add_argument("--tgs-minl", type=int, required=True, dest="minl", default=None, help="Min length of TGS reads.")
parser.add_argument("--protocol", type=str, required=True, dest="protocol", default=None, help="the error rate, error fraction and length distribution to the corresponding protocol will be used.")
args = parser.parse_args()


print("Simulating {} reads".format(args.protocol.upper()))

if args.protocol in ("ccs", "clr"):
    generate_PACBIO(args.pop_gen,
                    args.nread,
                    args.minl,
                    args.maxl,
                    args.outfa,
                    args.protocol)

if args.protocol in ("ont"):
    generate_ONT(args.pop_gen,
                 args.nread,
                 args.minl,
                 args.maxl,
                 args.outfa,
                 args.protocol)

print("Finished\n")

