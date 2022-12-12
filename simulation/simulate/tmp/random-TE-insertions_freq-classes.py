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


def parse_freqs(freqs):
    toret=[]
    if "," in freqs:
        toret=freqs.split(",")
    else:
        toret.append(freqs)
    return [float(i) for i in toret]
    
def create_insertiontuples(tefamcount,insertioncount,frequencies):
    toret=[]
    for i in range(1,tefamcount+1):
        for f in frequencies:
            for k in range(0,insertioncount):
                toret.append((i,f))
    random.shuffle(toret)
    return toret


def create_insertiontuplesrandfamfreq(tefamcount,insertioncount,frequencies):
    toret=[]
    for i in range(0,insertioncount):
        famnum=random.randint(0,tefamcount-1)
        freq=random.randint(0,len(frequencies)-1)
        toret.append((famnum,frequencies[freq]))
    random.shuffle(toret)
    return toret

def parse_subregion(subregion):
    if "-" not in subregion:
        raise Exception("Invalid subregion")
    substart,subend=subregion.split("-")
    return int(substart),int(subend)
    

parser.add_argument("--chasis", type=str, required=True, dest="ref_fasta", default=None, help="the chasis")
parser.add_argument("--te-seqs", type=str, required=True, dest="te_fasta", default=None, help="provide the TE consensus fasta file")
parser.add_argument("--pop-size", type=int, required=True, dest="pop_size", default=None, help="provide the number of haploid genomes")
parser.add_argument("--pop-freq", type=str, required=True, dest="pop_freq", default=None, help="provide a comma seperated list of population frequencies. i.e. \"0.1,0.3,0.5,0.6,1.0\"")
parser.add_argument("--insert-count", type=int, required=False, dest="insert_count", default=1, help="provide the number of insertions to create")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")
parser.add_argument("--min-distance", type=int, required=False, dest="min_distance", default=1000, help="Provide the minimum distance to next TE insertion")
parser.add_argument("--sub-region", type=str, required=False,dest="subregion",default=None,help="Restrict random insertion sites to a subregion of the chasis, e.g. 10000-20000")
args = parser.parse_args()



# load the genomes of the chasis and TEs
header,seq=load_chasis(args.ref_fasta)
fh=FastaReader.readFastaHash(args.te_fasta)


# compute the parameters
targetfrequencies=parse_freqs(args.pop_freq)

# subregion
substart=1
genomesize=len(seq)
subend=genomesize
if(args.subregion is not None):
    substart,subend=parse_subregion(args.subregion)
simulatesize=subend-substart
assert substart>0
assert simulatesize>0


tefamilycount=len(fh.keys())
insertions=create_insertiontuplesrandfamfreq(tefamilycount,args.insert_count,targetfrequencies)
popsize=args.pop_size

mindist=args.min_distance
maxdist=int(simulatesize/(len(insertions)+1))
if(maxdist<=mindist):
    raise ValueError("Genome of "+str(simulatesize)+ " is too small for "+str(len(insertions))+ " insertions with a min-distance between insertions of "+str(mindist))

tw = TEDefinitionWriter(args.output)

counter=substart
for famnum,freq in insertions:
    dist=random.randint(mindist,maxdist)
    counter+=dist
    inscount=int(freq*popsize)
    noninscount=popsize-inscount
    plusstrand=True
    if(random.random()<0.5):
        plusstrand=False
    
    tedefinitions=[TEInsertionSite(famnum,plusstrand) for i in range(0,inscount)]+[TEInsertionSite(0,True) for i in range(0,noninscount)]
    random.shuffle(tedefinitions) 
    tw.write(TEInsertionDefinition(str(counter),"tsd=0",tedefinitions))
    
    





