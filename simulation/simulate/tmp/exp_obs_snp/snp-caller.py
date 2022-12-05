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
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


from pileupIO import PileupCountReader



parser = argparse.ArgumentParser(description="""           
Description
-----------
Identifies SNPs, true positives (simulated and identified), false positives (not simulated but identified) and false negatives (simulated but not identified)
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output
----------
The SNPs, and some information to each SNP,
for example:
True	True	2R	201	T	36	1.0	0.472222222222	18:17:0:1
False	True	2R	250	C	34	1.0	0.941176470588	0:0:32:2
False	True	2R	398	T	36	0.0	0.944444444444	2:34:0:0

col 1: simulated SNP (True or False)
col 2: polymorphic SNP (True or False)
col 3: reference chromosome
col 4: position
col 5: reference character
col 6: coverage
col 7: strand bias; 0.5 = unbiased 1.0 = 100% forward strand, 0.0 = 100% reverse strand; 
col 8: frequency of reference allele
col 9: allele counts in the form A:T:C:G



Authors
-------
    Robert Kofler
""")


parser.add_argument("--pileup", type=str, required=True, dest="pileup", default=None, help="A pileup file")
parser.add_argument("--min-count", type=int, required=False, dest="mincount", default=2, help="the minimum count")
args = parser.parse_args()
mc=args.mincount

for p in PileupCountReader(args.pileup,0):
     
     positive=False
     polymorphic=False
     if(p.pos%100==1 and p.pos > 1):
          positive=True
     sample=p.samples[0]
     
     if sample.isPolymorphic(mc):
          polymorphic=True
     if positive or polymorphic:
          
     
          topr=[str(positive),str(polymorphic),p.chr,str(p.pos),p.refc]
          cov=sample.get_cov()
          crefc=sample.get_XXX(p.refc)
          cfwd=sample.get_forward()
          frefc="Na"
          ffwd="Na"
          if float(cov)>0.0:
               frefc = float(crefc)/float(cov)
               ffwd = float(cfwd)/float(cov)
          
          topr.append(str(cov))
          topr.append(str(ffwd))
          topr.append(str(frefc))
          topr.append("{0}:{1}:{2}:{3}".format(sample.get_sumA(),sample.get_sumT(),sample.get_sumC(),sample.get_sumG()) )
          print "\t".join(topr)

        
    
        
    

