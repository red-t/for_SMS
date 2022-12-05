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
from Multimeasure import *


def get_het(p):
     q=1.0-p
     h=1-p**2-q**2
     return h


def get_fst(ma1,mi1,ma2,mi2):
     if (ma1+mi1)<0.0000001:
          return 0.0
     if(ma2+mi2)<0.0000001:
          return 0.0
     p1=float(ma1)/float(ma1+mi1)
     p2=float(ma2)/float(ma2+mi2)
     
     h1=get_het(p1)
     h2=get_het(p2)
     hbet=(h1+h2)/2.0
     
     p12=(p1+p2)/2.0
     ht=get_het(p12)
     
     fst=(ht-hbet)/ht
     
     return fst


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
parser.add_argument("--min-count", type=int, required=False, dest="mincount", default=2, help="the minimum count")
parser.add_argument("--min-cov", type=int, required=False, dest="mincov", default=10, help="the minimum count")
parser.add_argument("--max-sb", type=float, required=False, dest="maxsb", default=0.5, help="the minimum count")
args = parser.parse_args()
mc=args.mincount






for s in StrandSyncReader(args.ssync):
     if not s.is_coveragesufficient(args.mincov):
          continue
     if s.is_strandbiased(args.maxsb):
          continue
     if not s.isPolymorphic(args.mincount):
          continue
     chr=s.chr
     pos=s.pos
     maj,min=s.get_maj_min()
     majc,minc=s.get_XXX(maj),s.get_XXX(min)
     fst=get_fst(majc[0],minc[0],majc[1],minc[1])
     
     topr=[chr,str(pos),"1","1.00","40.0"]
     topr.append("1:2={0}".format(fst))
     print "\t".join(topr)
     
     """2R	3722	1	1.000	21.0	1:2=0.05000000
2R	3729	1	1.000	38.0	1:2=0.03167421
2R	3736	1	1.000	48.0	1:2=0.02491103
2R	3738	1	1.000	53.0	1:2=0.03718200
2R	3739	1	1.000	56.0	1:2=0.02127660
     """


    
        
    

