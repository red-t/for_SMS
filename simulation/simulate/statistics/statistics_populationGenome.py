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




parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script computes statistics for the population genome""",formatter_class=argparse.RawDescriptionHelpFormatter,
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
parser.add_argument("--read-length", type=int, required=False, dest="read_length", default=None, help="Provide the read length to be simulates")
parser.add_argument("--reads", type=int, required=False, dest="reads", default=None, help="Provide the average physical coverage to simulate number of reads")
parser.add_argument("--coverage", type=int, required=False, dest="coverage", default=None, help="Provide the prefix for fastq files")
args = parser.parse_args()

reads=args.reads
readleng=args.read_length
coverage=args.coverage


totlen=0
totc=0
for h,s in FastaReader(args.pop_gen):
     rl=len(s)
     print "{0}\t{1}".format(h,rl)
     totlen+=rl
     totc+=1
avgenome=float(totlen)/float(totc)
print "Average genome length = {0}".format(avgenome)

if readleng is not None:
     if reads is not None:
          expectedcoverage=reads*2.0*readleng/avgenome
          print "Expected coverage with readleng {0} and read-count {1} = {2}".format(readleng,reads,expectedcoverage)
     
     if coverage is not None:
          reqreads=float(coverage*avgenome)/float(2.0*readleng)
          print "Required reads for a coverage {0} and using readlength {1} =  {2}".format(coverage,readleng,reqreads)
          




