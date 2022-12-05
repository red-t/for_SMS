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
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


from pileupIO import StrandSyncReader


def get_sync(ss):
     toret="{0}:{1}:{2}:{3}:0:0".format(ss.get_sumA(),ss.get_sumT(),ss.get_sumC(),ss.get_sumG())
     return toret




parser = argparse.ArgumentParser(description="""           
Description
-----------
Identifies SNPs, true positives (simulated and identified), false positives (not simulated but identified) and false negatives (simulated but not identified)
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")


parser.add_argument("--ssync", type=str, required=True, dest="ssync", default=None, help="A pileup file")
parser.add_argument("--min-cov", type=int, required=False, dest="mincov", default=10, help="the minimum coverage")
parser.add_argument("--max-sb", type=float, required=False, dest="maxsb", default=0.5, help="the maximum strand bias")
args = parser.parse_args()



for s in StrandSyncReader(args.ssync):
     if not s.is_coveragesufficient(args.mincov):
          continue
     if s.is_strandbiased(args.maxsb):
          continue

     chr=s.chr
     pos=s.pos
     refc=s.refc
     topr=[chr,str(pos),refc]
     for i in s.samples:
          topr.append(get_sync(i))
     print "\t".join(topr)
          

     
     """
     2R	2	A	2:0:0:0:0:0	2:0:0:0:0:0
     2R	3	A	4:0:0:0:0:0	4:0:0:0:0:0
     2R	4	G	0:0:0:5:0:0	0:0:0:5:0:0
     2R	5	T	0:5:0:0:0:0	0:5:0:0:0:0
     """


    
        
    

