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
import Mutator
import CoverageGenerator
import ReadLengthDistribution



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
parser.add_argument("--read-length", type=int, required=False, dest="read_length", default=None, help="the mean read length")
parser.add_argument("--tgs-alpha", type=float, required=False, dest="tgs_alpha", default=None, help="alpha of gamma model applied for TGS length generation.")
parser.add_argument("--tgs-loc", type=float, required=False, dest="tgs_loc", default=None, help="loc of gamma model applied for TGS length generation.")
parser.add_argument("--tgs-beta", type=float, required=False, dest="tgs_beta", default=None, help="beta of gamma model applied for TGS length generation.")
parser.add_argument("--rld-file", type=str, required=False, dest="rldfile", default=None, help="read length distribution file; will override --read-length and --std-dev")
parser.add_argument("--error-rate", type=float, required=False, dest="error_rate", default=0.0, help="the error rate of the reads; indels")
parser.add_argument("--error-fraction", type=str, required=False, dest="err_frac", help="the fraction of each type of errors, mis:del:ins")
parser.add_argument("--reads", type=int, required=True, dest="reads", default=None, help="the total number of reads")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="output - a fasta file")
parser.add_argument("--tgs-maxl", type=int, required=True, dest="tgs_maxl", default=None, help="Max length of TGS reads")
parser.add_argument("--tgs-minl", type=int, required=True, dest="tgs_minl", default=None, help="Min length of TGS reads.")



args = parser.parse_args()

print "Simulating PacBio Pool-Seq reads"

print "Reading the length of the population genome"
pgld=fastaIO.SequenceUtility.get_length_list(args.pop_gen)

print "Getting mutator for PacBio reads with error rate {0} and error fraction {1}".format(args.error_rate,args.err_frac)
mutator=Mutator.PacBioMutator(args.error_rate,args.err_frac) # get a suitable mutator; suitability depends on the error rate

readnumbergenerator=CoverageGenerator.RandomReads(args.reads,pgld)

rldfactory=ReadLengthDistribution.get_rld_factory(args.read_length, args.tgs_alpha, args.tgs_loc, args.tgs_beta, args.rldfile)

# get single-end writer
fawriter=fastaIO.FastaWriter(args.fasta,60)

counter=0
readcount=1
f_err = open("ErrorPerRead.txt", "a")
for header,seq in fastaIO.FastaReader(args.pop_gen):
    targetreads=readnumbergenerator.get_reads(counter)
    print "Generating {0} reads for haploid genome {1}".format(targetreads,header)
    counter+=1
    
    sl=len(seq)
    for i in range(0,targetreads):
        readlen=rldfactory.next()
        if readlen < args.tgs_minl:
            readlen = random.randint(args.tgs_minl, args.tgs_maxl)
        if readlen > args.tgs_maxl:
            readlen = random.randint(args.tgs_minl, args.tgs_maxl)
            
        firstposition=random.randint(0,sl-readlen)
        read1=seq[firstposition:firstposition+readlen]
        if random.random()<0.5:
               read1=fastaIO.SequenceUtility.rc(read1)
        
        read1, errs=mutator.mutateseq(read1)
        
        f_err.write('{0};{1}:{2}\t{3}\t{4}\t{5}\n'.format(readcount,header,firstposition,errs[0],errs[1],errs[2]))
        h="{0};{1}:{2}".format(readcount,header,firstposition)
        fawriter.write(h,read1)
        readcount+=1
    

f_err.close()
fawriter.close()
print "Finished"

