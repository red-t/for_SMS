from NGS_utils import generate_PE
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--pg", type=str, required=True, dest="pop_gen", default=None, help="the population genome - a fasta file")
parser.add_argument("--read-length", type=int, required=True, dest="read_length", default=None, help="the read length")
parser.add_argument("--inner-distance", type=int, required=True, dest="inner_distance", default=None, help="the inner distance between the paired end reads")
parser.add_argument("--std-dev", type=int, required=True, dest="std_dev", default=None, help="the standard deviation of the inner distance")
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="the error rate of the reads")
parser.add_argument("--reads", type=int, required=True, dest="reads", default=None, help="the total number of paired-end reads")
parser.add_argument("--fraction-chimera", type=float, required=False, dest="chimera", default=0.0, help="fraction of chimeric paired ends; i.e. paired-ends from unrelated positions")
parser.add_argument("--fastq1", type=str, required=True, dest="fastq1", default=None, help="output fastq file - first read")
parser.add_argument("--fastq2", type=str, required=True, dest="fastq2", default=None, help="output fastq file - second read")
args = parser.parse_args()

print("Simulating Illumina paired-end reads for Pool-Seq")
print("Getting mutator for Illumina reads with error rate {0}; Base substitutions with Poisson distributed errors".format(args.error_rate))

generate_PE(args.pop_gen,
            args.error_rate,
            args.reads,
            args.inner_distance,
            args.std_dev,
            args.read_length,
            args.fastq1,
            args.fastq2,
            args.chimera)

print("Finished")