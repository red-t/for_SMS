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
    This script simulates single-end reads from the population genome""",formatter_class=argparse.RawDescriptionHelpFormatter,
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
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="the error rate of the reads")
parser.add_argument("--reads", type=int, required=True, dest="reads", default=None, help="the total number of single-end reads")
parser.add_argument("--fastq", type=str, required=True, dest="fastq", default=None, help="output fastq file")



args = parser.parse_args()

print "Simulating Illumina single-end reads for Pool-Seq"

print "Reading the length of the population genome"
pgld=fastaIO.SequenceUtility.get_length_list(args.pop_gen)

print "Getting mutator for Illumina reads with error rate {0}; Base substitutions with Poisson distributed errors".format(args.error_rate)
mutator=Mutator.PoisonSeqMutator(args.error_rate) # get a suitable mutator; suitability depends on the error rate

readnumbergenerator=CoverageGenerator.RandomReads(args.reads,pgld)

readleng=args.read_length

# get single-end writer
fqwriter=fastqIO.FastqWriter(args.fastq)

counter=0
readcount=1
for header,seq in fastaIO.FastaReader(args.pop_gen):
    targetreads=readnumbergenerator.get_reads(counter)
    print "Generating {0} reads for haploid genome {1}".format(targetreads,header)
    counter+=1
    
    sl=len(seq)
    for i in range(0,targetreads):
        
        firstposition=random.randint(0,sl-readleng)
        read1=seq[firstposition:firstposition+readleng]
        if random.random()<0.5:
               read1=fastaIO.SequenceUtility.rc(read1)
        read1=mutator.mutateseq(read1)
        
        h="{0};{1}:{2}".format(readcount,header,firstposition)
        fqwriter.write(h,read1)
        readcount+=1
    

fqwriter.close()
print "Finished"

