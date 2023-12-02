import argparse
from define_pgdf_utils import getGermOrder, defineHeader, defineBody


def getArgs():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument("--chassis", type=str, dest="templateFa", default=None, help="the template, i.e. the sequence into which TEs will be inserted; a fasta file")
    parser.add_argument("--te-seqs", type=str, dest="teFa", default=None, help="TE CSS sequences in a fasta file")
    parser.add_argument("--N", type=int, dest="numTotalGenome", default=None, help="the number of haploid genomes")
    parser.add_argument("--divergence-rate", type=float, dest="insDivRate", help="TE sequence divergence rate")
    parser.add_argument("--germline-count", type=int, dest="numGerm", default=1, help="germline insertions mumber in the population genome")
    parser.add_argument("--somatic-count", type=int, dest="numSoma", default=1, help="somatic insertions mumber in the population genome")
    parser.add_argument("--output", type=str, dest="output", default=None, help="the output file")
    parser.add_argument("--min-distance", type=int, dest="minDist", default=1000, help="minimum distance between TE insertions")
    parser.add_argument("--species", type=str, dest="species", default="human", help="the transposon library to the corresponding species will be used")
    parser.add_argument("--truncProb", type=float, dest="truncProb", help="Probability of constructing a truncation on an inserted sequence")
    parser.add_argument("--nestProb", type=float, dest="nestProb", help="Probability of constructing a nested insertion on an inserted sequence")
    parser.add_argument("--mode", type=int, dest="mode", help="Generate data for: (1)model training or (2)benchmarking")
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = getArgs()

    if (args.numGerm + args.numSoma) > 0:        
        germOrderSet = getGermOrder(args.numGerm, args.numSoma)

        # define header of the pgd-file (population genome definition)
        insIdToExpressDict, chrom, chromLen, insIdList = defineHeader(args, germOrderSet)
        
        # define body of the pgd-file (population genome definition)
        minDist = args.minDist
        try:
            maxDist = int(chromLen / (args.numGerm + args.numSoma))
        except:
            maxDist = int(chromLen / (args.numGerm + args.numSoma + 1))
        
        if(maxDist <= minDist):
            raise ValueError("Genome of "+str(chromLen)+ " is too small for "+str(args.numGerm + args.numSoma)+ " insertions with a min-distance between insertions of "+str(minDist))

        defineBody(insIdToExpressDict, chrom, minDist, maxDist, germOrderSet, insIdList, args)