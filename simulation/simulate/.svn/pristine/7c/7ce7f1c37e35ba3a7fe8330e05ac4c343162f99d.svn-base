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

from fastaIO import FastaReader
from FastqPairWriter import *
import SeqLoader
import Mutator
from TargetCoverage import get_uniform_readnumbergenerator




parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script simulate paired-end reads from the population genomes based on various parameters""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    (1) python version 3.4.3  

Authors
-------
    Ram Vinay Pandey
    Robert Kofler
    Christian Schloetterer    
""")


parser.add_argument("--population-genome", type=str, required=True, dest="pop_gen", default=None, help="Provide the population genome fasta file created by 03_create_population_genome.py")
parser.add_argument("--read-length", type=int, required=True, dest="read_length", default=None, help="Provide the read length to be simulates")
parser.add_argument("--inner-distance", type=int, required=True, dest="inner_distance", default=None, help="Provide the inner distance between forward and reverse read")
parser.add_argument("--std-dev", type=int, required=True, dest="std_dev", default=None, help="Provide the standard deviation inner distance between read two reads")
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="Provide the sequencing error rate to be simulated")
parser.add_argument("--coverage", type=str, required=False, dest="coverage", default=None, help="either 2.0pcpg average physical coverage per genome, 2.0cpg average coverage per genome, 1000rpg reads per genome")
parser.add_argument("--fastq-prefix", type=str, required=True, dest="fastq_prefix", default=None, help="Provide the prefix for fastq files")
#parser.add_argument("--individual-sequencing", type=bool,required=False, dest="individual_seq",default=False,help="Switch to sequencing individuals separately; default is the Pool-seq approach; default=False")
#parser.add_argument("--fraction-chimera", type=float, required=False, dest="chimera", default=0.0, help="Fraction of chimera")

args = parser.parse_args()

print "Generating paired end reads"
insertsize=args.inner_distance
stddev=args.std_dev
print "Loading the length of the population genome"
pgll=SeqLoader.getPopGenomeLengthList(args.pop_gen)

print "Getting uniform read number generator"
readnumbergenerator=get_uniform_readnumbergenerator(args.coverage ,args.inner_distance ,args.read_length, pgll)


readleng=args.read_length
if(args.error_rate>=1.0):
     raise Exception("Error rate is to high")
mutator=Mutator.getSeqMutator(args.error_rate) # get a suitable mutator; suitability depends on the error rate
fw=FastqPoolWriter(args.fastq_prefix)

c=0
for header,seq in FastaReader(args.pop_gen):
    sl=len(seq)
    reads=readnumbergenerator.get_reads(c)
    c+=1
    avoffset=float(sl)/float(reads)
    pos=0.0
    print "Lenght of {0} = {1}".format(header,sl)
    print "Number of reads {0}".format(reads)
    for i in range(0,reads):
        firstposition=int(pos)
        innerdistance=int(random.gauss(insertsize,stddev))
        while(innerdistance<1):
          innerdistance=int(random.gauss(insertsize,stddev))
        secondposition=firstposition+readleng+innerdistance
        
        if(secondposition+readleng>=sl):
          continue
        read1=seq[firstposition:firstposition+readleng]
        read2=SeqLoader.rc(seq[secondposition:secondposition+readleng])
        if(len(read1)!=readleng):
          raise Excpetion("Invalid read {0} at position {1}".format(read1,firstposition))
        if(len(read2)!=readleng):
          raise Excpetion("Invalid read {0} at position {1}".format(read2,secondposition))        
        read1=mutator.mutateseq(read1)
        read2=mutator.mutateseq(read2)
        fw.write(read1,read2,i+1)
        pos+=avoffset #float level computation
fw.close()
print "Finished"

