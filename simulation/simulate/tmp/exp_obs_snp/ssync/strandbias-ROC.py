#!/usr/bin/env python

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

from Multimeasure import *


totol=1000.0

parser = argparse.ArgumentParser(description="""           

Identifies SNPs, true positives (simulated and identified), false positives (not simulated but identified) and false negatives (simulated but not identified)
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")
parser.add_argument("--sb", type=str, required=True, dest="sbfile", default=None, help="A strandbias file")
parser.add_argument("--points", type=int, required=True, dest="points", default=None, help="the points")
args = parser.parse_args()



vals=[]
for m in MultimeasureReader(args.sbfile):
     if not m.issnp:
         continue
     vals.append(m)
     
count=len(vals)
steps=count/args.points

vals=sorted(vals,key=lambda m: -m.strandbias)

nonol=float(len(vals)-totol)

count=0
notoutliercount=0
outliercount=0
for v in vals:
     count+=1
     if v.isoutlier:
          outliercount+=1
     else:
          notoutliercount+=1
          
     if count>1 and count % steps==0:

          print "{0}\t{1}\t{2}".format(float(outliercount)/totol,float(notoutliercount/nonol),v.strandbias)
          
     


     

