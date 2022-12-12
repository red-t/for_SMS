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


from pileupIO import PileupCountReader,StrandSyncWriter



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


parser.add_argument("--pileup", type=str, required=True, dest="pileup", default=None, help="A pileup file")
parser.add_argument("--min-qual", type=int, required=False, dest="minqual", default=0, help="the minimum quality")
args = parser.parse_args()



sw=StrandSyncWriter("/dev/stdout") 



for p in PileupCountReader(args.pileup,args.minqual):
     sw.write(p)
sw.close()

    
        
    

