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


parser = argparse.ArgumentParser(description="""
Description
-----------
This script generates an empty TE landscape for manual editing
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")

from fastaIO import FastaReader, FastaWriter, SequenceUtility
from PopGenomeDefinitionIO import *
from TESequenceBuilder import SequenceContainer
from TEInsert import *



 
    
parser.add_argument("--chassis", type=str, required=False, dest="ref_fasta", default=None, help="the chassis, i.e. the sequence into which TEs will be inserted; a fasta file")
parser.add_argument("--te-seqs", type=str, required=False, dest="te_fasta", default=None, help="TE sequences in a fasta file")
parser.add_argument("--pgd", type=str, required=True, dest="pgd_definition", default=None, help="the definition of the population genome")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file; will be multi-fasta file")
parser.add_argument("--ins-seq", type=str, required=False, dest="ins_seq", default=None, help="the output file of insertion sequence; will be multi-fasta file")
parser.add_argument("--sub_idx", type=int, required=True, dest="sub_idx", default=None, help="the index of the sub-population genome")
parser.add_argument("--sub_size", type=int, required=True, dest="sub_size", default=None, help="the size of the sub-population genome")
args = parser.parse_args()

# read TE sequences from file; if provided
tetuples=[]
if args.te_fasta is not None:
     tmp=FastaReader.readAllTuples(args.te_fasta)
     tetuples=[t[1] for t in tmp]
     print "Loading TE sequences; Found {0} in file {1}".format(len(tetuples),args.te_fasta)
sc=SequenceContainer(tetuples)

# read the PGD; must be provided
print "Loading population genome defintion"
pgdr=PopGenDefinitionReader(args.pgd_definition,sc) # construct sequence(s) which will be inserted into chasis, mainly change pgdr.tuples, sc.__sc
tedeftuples=pgdr.read_transposed() # extract insertions that will be inserted to each haploid genome, [[(pos1, insid1),(pos2, insid2), ...], [(pos1, insid1)], ...]
print "Found {0} TE defintions".format(sc.get_count_definitions())
print "Will simulate {0} TE insertion sites within a population having {1} haploid genomes".format(pgdr.insertions, pgdr.popsize)

# load chasis from the file; if provided otherwise from the PGD; not both though
chasis=""
if args.ref_fasta is not None:
     if pgdr.get_chasis() !="":
          raise Exception("Two chasis were provided (fasta file and pop genome definition); invalid, provide only one")
     print "Loading chasis from file "+args.ref_fasta
     crap,chasis=SequenceUtility.load_chasis(args.ref_fasta) # seq_id, sequence
else:
     print "No chasis file found; Will use chasis from the population definition file"
     crap="TOY"
     chasis=pgdr.get_chasis()
if chasis=="":
     raise Exception("No chasis was provided, neither in a fastq file nor in the population genome definition file")
print "Will proceed with chasis having a size of {0} nt".format(len(chasis))






counter=1 + args.sub_idx*args.sub_size
fw=FastaWriter(args.output,60)
print "Start writing population genome"
for tmp in tedeftuples: # [(pos1, insid1),(pos2, insid2), ...]
     # translate into sequences
     seqidtup = [(t[0],str(crap)+"_hg"+str(counter)+";"+t[1]) for t in tmp]
     seqtup=[(t[0],sc.getTESequence(t[1])) for t in tmp] # [(pos1, TESequence),(pos2, TESequence), ...]
     seq_with_te=SeqInserter.insertSequences(chasis,seqtup, args.ins_seq, seqidtup)
     # TODO seqinserter must check if position is smaller than sequence length

     fw.write(str(crap)+"_hg"+str(counter),seq_with_te)
     counter+=1
print "Done; Wrote {0} genomes to file {1}".format(args.sub_size,args.output)
fw.close()


# if args.ins_seq:
#      seqid=list(set(seqid))
#      fout=open(args.ins_seq, "a")
#      for id in seqid:
#           fout.write("{0}\t{1}\n".format(str(id), sc.getTESequence(id).sequence))

#      fout.close()
