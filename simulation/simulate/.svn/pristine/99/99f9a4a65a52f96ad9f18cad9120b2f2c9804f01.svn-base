#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import collections


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
from TEDefinition import TEDefinitionReader,TEInsertionSite
from TEHierarchy import *

parser = argparse.ArgumentParser(description="""           
Description
-----------
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Ram Vinay Pandey
    Christian Schloetterer
 
""")

def loadnum2famtrans(teseqar,hier):
    n2f={}
    for i in range(1,len(teseqar)+1):
        fam=teseqar[i-1]
        fam=hier.getFam(fam)
        n2f[i]=fam
    return n2f


def load_teseqsnamearray(tefile):
    toret=[]
    fr=FastaReader(tefile)
    for header,seq in fr:
        toret.append(header)
    return toret        
    

parser.add_argument("--te-sites", type=str, required=True, dest="te_sites", default=None, help="A file with the TE insertion sites")
parser.add_argument("--te-seqs", type=str, required=True, dest="te_seqs", default=None, help="A file with the TE sequences")
parser.add_argument("--hier", type=str, required=False, dest="hier", default="", help="A TE hierarchy; if provided output will be for families")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="The output file")
args = parser.parse_args()

# number 2 family name translator
hier= TEHierarchyDefault()
if args.hier !="":
    hier=loadtehier(args.hier)
teseqarray=load_teseqsnamearray(args.te_seqs)
num2fam=loadnum2famtrans(teseqarray,hier)

tesites= TEDefinitionReader.readall(args.te_sites)
ofh=open(args.output,"w")

for s in tesites:
    p=s.position
    fcounter=collections.defaultdict(lambda:0)
    sites=0
    for tes in s.tesites:
        fcounter[tes.getKey()]+=1
        sites+=1
    ftem=fcounter.items()
    ftem=sorted(ftem,key=lambda k:-k[1])
    
    topr=[]
    topr.append(str(p))
    for tekey,count in ftem:
        teids=TEInsertionSite.parse_key(tekey)
        teid=teids.teid
        if(teid==0):
            continue
        fam=num2fam[teid]
        freq=float(count)/float(sites)
        topr.append(fam)
        if(teids.plusstrand):
            topr.append("+")
        else:
            topr.append("-")
        topr.append(str(freq))
    ofh.write("\t".join(topr)+"\n")
        
    
        
    

