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
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Ram Vinay Pandey
    Christian Schloetterer
 
""")


parser.add_argument("--pileup", type=str, required=True, dest="pileup", default=None, help="A pileup file")

args = parser.parse_args()

for p in PileupCountReader(args.pileup,20):
     pos=p.pos
     if(pos%100!=1):
          continue
     topr=[p.chr,str(pos),p.refc]
     s=p.samples[0]
     cov=s.get_cov()
     crefc=s.get_XXX(p.refc)
     cfwd=s.get_forward()
     frefc=float(crefc)/float(cov)
     ffwd=float(cfwd)/float(cov)
     
     topr.append(str(cov))
     topr.append(str(ffwd))
     topr.append(str(frefc))
     print "\t".join(topr)

        
    
        
    

