#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import collections
import math


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


def update_list(toplist,top,key,me):
     if(len(toplist)<top):
          toplist.append((me,key))
     else:
          toplist=sorted(toplist)
          if(toplist[0][0] < me):
               toplist[0]=(me,key)
     assert len(toplist)<=top
     return toplist



parser = argparse.ArgumentParser(description="""           
Possible strandbiases
a.) maxfraction, the maximum strandbias in any population; not cosidering allele frequency differences


-----------
Identifies SNPs, true positives (simulated and identified), false positives (not simulated but identified) and false negatives (simulated but not identified)
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")
parser.add_argument("--mm", type=str, required=True, dest="mm", default=None, help="A pileup file")
parser.add_argument("--top", type=int, required=True, dest="top", default=None, help="the top-n outlier to mark")
args = parser.parse_args()




toplist=[]
mr=MultimeasureReader(args.mm)
for mm in mr:
     key="{0}:{1}".format(mm.chr,mm.pos)
     me=mm.measure
     toplist=update_list(toplist,args.top,key,me)
mr.close()
keyset=set([i[1] for i in toplist])

mr=MultimeasureReader(args.mm)
w=MultimeasureWriter("/dev/stdout")
for mm in mr:
     key="{0}:{1}".format(mm.chr,mm.pos)
     if key in keyset:
          mm.isoutlier=True
     else:
          mm.isoutlier=False
     w.write(mm)
mr.close()
w.close()

     


