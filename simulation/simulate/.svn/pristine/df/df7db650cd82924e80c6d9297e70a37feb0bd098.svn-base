#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import inspect


 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)

import fastaIO
import fastqIO
import Mutator
import CoverageGenerator



parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script simulate paired-end reads from the population genomes based on various parameters""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+ 

Authors
-------
    Robert Kofler  
""")


parser.add_argument("--pg", type=str, required=True, dest="pop_gen", default=None, help="the population genome - a fasta file")
parser.add_argument("--read-length", type=int, required=True, dest="read_length", default=None, help="the read length")
parser.add_argument("--inner-distance", type=int, required=True, dest="inner_distance", default=None, help="the inner distance between paired-end reads")
parser.add_argument("--std-dev", type=int, required=True, dest="std_dev", default=None, help="the standard deviation of the inner distance")
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="the error rate of the reads")
parser.add_argument("--reads", type=int, required=False, dest="reads", default=None, help="the total number of paired-end reads per individual (diploid or haploid)")
parser.add_argument("--fraction-chimera", type=float, required=False, dest="chimera", default=0.0, help="fraction of chimeric paired ends; i.e. pairs from unrelated positions")
parser.add_argument("--fastq-prefix", type=str, required=True, dest="fastq_prefix", default=None, help="output - prefix of the fastq files")
parser.add_argument("--haploid",action="store_true", dest="haploid", help="Switch to haploid genomes; if not provided diploid genomes are used; default=False")

args = parser.parse_args()

print "Simulating Illumina paired-end reads for sequencing individuals"


print "Getting mutator for Illumina reads with error rate {0}; Base substitutions with Poisson distributed errors".format(args.error_rate)
mutator=Mutator.PoisonSeqMutator(args.error_rate) # get a suitable mutator; suitability depends on the error rate



insertsize=args.inner_distance
stddev=args.std_dev
readleng=args.read_length

# get paired-end writer
fqwriter=fastqIO.FastqPEBatchWriter(args.fastq_prefix,args.haploid)

reads=args.reads
if(not args.haploid):
     reads=int(reads/2.0)

counter=0
readcount=1
for header,seq in fastaIO.FastaReader(args.pop_gen):
    targetreads=reads
    print "Generating {0} reads for haploid genome {1}".format(targetreads,header)
    counter+=1
    
    sl=len(seq)
    for i in range(0,targetreads):
        firstposition=None
        secondposition=None
        ischimera=random.random()
        if(ischimera<args.chimera): # read is a chimera
            firstposition=random.randint(0,sl-readleng)
            secondposition=random.randint(0,sl-readleng)
        else:                   # read is normal paire-end
            innerdistance=int(random.gauss(insertsize,stddev))
            while(innerdistance<1): # must be larger than 0
               innerdistance=int(random.gauss(insertsize,stddev))
            outerdistance=2*readleng+innerdistance
            firstposition=random.randint(0,sl-outerdistance)
            secondposition=firstposition+readleng+innerdistance
        read1=seq[firstposition:firstposition+readleng]
        read2=fastaIO.SequenceUtility.rc(seq[secondposition:secondposition+readleng])
        read1=mutator.mutateseq(read1)
        read2=mutator.mutateseq(read2)
        
        h="{0};{1}:{2}-{3}".format(readcount,header,firstposition,secondposition)
        fqwriter.write(h,read1,read2,counter)
        readcount+=1
    

fqwriter.close()
print "Finished"

