#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random


 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


from SeqLoader import *
from fastaIO import FastaReader,FastaWriter
from TEDefinition import TEDefinitionWriter,TEInsertionSite,TEInsertionDefinition

parser = argparse.ArgumentParser(description="""           
Description
-----------
This script generates random insertions of TE insertions
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Ram Vinay Pandey
    Robert Kofler
    Christian Schloetterer
 
""")

    

parser.add_argument("--chasis", type=str, required=True, dest="ref_fasta", default=None, help="the chasis")
parser.add_argument("--te-seqs", type=str, required=True, dest="te_fasta", default=None, help="provide the TE consensus fasta file")
parser.add_argument("--pop-size", type=int, required=True, dest="pop_size", default=None, help="provide the number of haploid genomes")
parser.add_argument("--freq-start", type=float, required=True, dest="freq_start", default=None, help="start of the frequency range")
parser.add_argument("--freq-end", type=float, required=True, dest="freq_end", default=None, help="end of the frequency range")
parser.add_argument("--insert-count", type=int, required=False, dest="insert_count", default=1, help="provide the total number of insertions to create")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")
parser.add_argument("--min-distance", type=int, required=False, dest="min_distance", default=1000, help="Provide the minimum distance to next TE insertion")
args = parser.parse_args()

# load the genomes of the chasis and TEs
header,seq=load_chasis(args.ref_fasta)
fh=FastaReader.readFastaHash(args.te_fasta)

# compute the parameters
genomesize=len(seq)
tefamilycount=len(fh.keys())
insertions=args.insert_count
popsize=args.pop_size

mindist=args.min_distance
maxdist=int(genomesize/(insertions+1))
if(maxdist<=mindist):
    raise ValueError("Genome of "+str(genomesize)+ " is too small for "+str(insertions)+ " insertions with a min-distance between insertions of "+str(mindist))

tw = TEDefinitionWriter(args.output)

counter=1
for i in range(0,insertions):
    dist=random.randint(mindist,maxdist)
    counter+=dist
    famnum=random.randint(1,len(fh))
    freq=random.uniform(args.freq_start,args.freq_end)
    plusstrand=True
    if(random.random()<0.5):
        plusstrand=False
    inscount=int(freq*popsize)
    noninscount=popsize-inscount
    
    # get the TE definition
    tedefinitions=[TEInsertionSite(famnum,plusstrand) for i in range(0,inscount)]+[TEInsertionSite(0,True) for i in range(0,noninscount)]
    random.shuffle(tedefinitions) 
    tw.write(TEInsertionDefinition(str(counter),"tsd=0",tedefinitions))
    
    





