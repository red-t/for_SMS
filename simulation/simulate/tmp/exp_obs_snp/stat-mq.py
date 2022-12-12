#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections
import fileinput




parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--win", type=int, required=False, dest="win", default=5, help="The window size of the mapping qualities")
args = parser.parse_args()
win=args.win

count_reads=0
count_unmapped=0
count_mapped=0
count_pp=0
mqh=collections.defaultdict(lambda:0)
mqwinh=collections.defaultdict(lambda:0)
for line in args.sam:
     """
     0  1        2      3       4         5     6        7       8
     1	99	2R	1	0	100M	=	201	300	
     """
     a=line.rstrip("\n").split("\t")
     flag=int(a[1])
     count_reads+=1
     if flag & 0x004 > 0:
          count_unmapped+=1
          continue
     count_mapped+=1
     mq=int(a[4])
     mqh[mq]+=1
     key=int(mq/win)
     mqwinh[key]+=1
     is_pp=flag & 0x0002 > 0
     mate_unmapped= flag & 0x0008 > 0     
     if is_pp and not mate_unmapped:
          count_pp+=1
          

print "Count\treads\t{0}".format(count_reads)
print "Count\tmapped\t{0}".format(count_mapped)
print "Count\tunmapped\t{0}".format(count_unmapped)
print "Count\tproperpair\t{0}".format(count_pp)


for i in sorted(mqwinh.keys()):
     print "MQwindows\t{0}\t{1}".format(i*win,mqwinh[i])
          
for i in sorted(mqh.keys()):
     print "MQdetailed\t{0}\t{1}".format(i,mqh[i])
     
     
