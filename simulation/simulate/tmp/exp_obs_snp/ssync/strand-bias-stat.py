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
parser.add_argument("--step", type=float, required=True, dest="step", default=None, help="the step size")
parser.add_argument("--min-coverage", type=int, required=False, dest="mincov", default=10, help="the step size")


args = parser.parse_args()

sbh=collections.defaultdict(lambda:collections.defaultdict(lambda:0 ))
sbm=collections.defaultdict(lambda:0)
totc=0

for ss in StrandSyncReader(args.ssync):
     if not ss.is_coveragesufficient(args.mincov):
          continue
     totc+=1
     maxsbk=0
     for i,sa in enumerate(ss.samples):
          sb=sa.get_strandbias()
          key=int(sb/args.step)
          if key>maxsbk:
               maxsbk=key
          sbh[key][i]+=1
     sbm[maxsbk]+=1


print "tot\t{0}".format(totc)
for cla in sorted(sbh.keys()):
     freqcl=cla*args.step
     topr=["sep",str(freqcl)]
     for s,count in sbh[cla].items():
          topr.append(str(count))
     print "\t".join(topr)


retro=0
for cla in reversed(sorted(sbm.keys())): 
     freqcl=cla*args.step
     count=sbm[cla]
     retro+=count
     fraction=(float(retro)/float(totc))
     topr=["max",str(freqcl),str(count),str(fraction)]

     print "\t".join(topr)
     
     
     
          

    
        
    

