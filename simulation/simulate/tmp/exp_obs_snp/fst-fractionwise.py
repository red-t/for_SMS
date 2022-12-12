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





parser = argparse.ArgumentParser(description="""           
Description
-----------
Input: requires
A FST file created from a sync file with 4-columns
     2R	3681	1	1.000	10.0	1:2=0.05263158	1:3=0.04347826	1:4=0.04347826	2:3=0.03336046	2:4=0.03336046	3:4=0.00000000
     2R	3685	1	1.000	15.0	1:2=0.00000000	1:3=0.03448276	1:4=0.03448276	2:3=0.03448276	2:4=0.03448276	3:4=0.00000000
     2R	3688	1	1.000	18.0	1:2=0.00000000	1:3=0.02857143	1:4=0.02857143	2:3=0.02857143	2:4=0.02857143	3:4=0.00000000
     2R	3691	1	1.000	19.0	1:2=0.03952883	1:3=0.01411503	1:4=0.01411503	2:3=0.02702703	2:4=0.02702703	3:4=0.00000000
     2R	3692	1	1.000	20.0	1:2=0.04892547	1:3=0.05130316	1:4=0.05130316	2:3=0.00000000	2:4=0.00000000	3:4=0.00000000

""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")


def parse_comparisions(compstr):
     if ":" not in compstr:
          raise Exception("At least one comparision needs to be provided")
     
     toparse=[]
     if "," in compstr:
          toparse=compstr.split(",")
     else:
          toparse=[compstr]
     return toparse

def parse_steps(steps,filesize):
     toret=[]
     if "," in steps:
          toret=steps.split(",")
     else:
          toret=[steps]
     
     floats=[float(t) for t in toret]
     ints=[int(filesize*t) for t in floats]

     return set(ints)



def parseline(line):
     """
     2R	3681	1	1.000	10.0	1:2=0.05263158	1:3=0.04347826	1:4=0.04347826	2:3=0.03336046	2:4=0.03336046	3:4=0.00000000
     2R	3685	1	1.000	15.0	1:2=0.00000000	1:3=0.03448276	1:4=0.03448276	2:3=0.03448276	2:4=0.03448276	3:4=0.00000000
     2R	3688	1	1.000	18.0	1:2=0.00000000	1:3=0.02857143	1:4=0.02857143	2:3=0.02857143	2:4=0.02857143	3:4=0.00000000
     2R	3691	1	1.000	19.0	1:2=0.03952883	1:3=0.01411503	1:4=0.01411503	2:3=0.02702703	2:4=0.02702703	3:4=0.00000000
     2R	3692	1	1.000	20.0	1:2=0.04892547	1:3=0.05130316	1:4=0.05130316	2:3=0.00000000	2:4=0.00000000	3:4=0.00000000
     """
     a=line.split("\t")
     chr=a.pop(0)
     pos=int(a.pop(0))
     t=a.pop(0)
     t=a.pop(0)
     t=a.pop(0)
     fsth={}
     for t in a:
          key,val=t.split("=")
          fsth[key]=float(val)
     return pos,fsth
     

          


parser.add_argument("--fst", type=str, required=True, dest="fst", default=None, help="A pileup file")
parser.add_argument("--fractions", type=str, required=False, dest="fractions", default="0.001,0.01,0.1", help="The numbers of windows")
parser.add_argument("--comparisions",type=str,required=False, dest="comparisions",default="1:2",help="the comparisions to perform")
args = parser.parse_args()

comp=parse_comparisions(args.comparisions)



positivesize=0
negativesize=0


tp=[]
fp=[]

for l in open(args.fst):
     a=l.rstrip("\n").split("\t")
     """
     0   1      2        3         4         5              6              7              8              9                   10
     2R	3681	1	1.000	10.0	1:2=0.05263158	1:3=0.04347826	1:4=0.04347826	2:3=0.03336046	2:4=0.03336046	3:4=0.00000000
     2R	3685	1	1.000	15.0	1:2=0.00000000	1:3=0.03448276	1:4=0.03448276	2:3=0.03448276	2:4=0.03448276	3:4=0.00000000
     2R	3688	1	1.000	18.0	1:2=0.00000000	1:3=0.02857143	1:4=0.02857143	2:3=0.02857143	2:4=0.02857143	3:4=0.00000000
     2R	3691	1	1.000	19.0	1:2=0.03952883	1:3=0.01411503	1:4=0.01411503	2:3=0.02702703	2:4=0.02702703	3:4=0.00000000
     2R	3692	1	1.000	20.0	1:2=0.04892547	1:3=0.05130316	1:4=0.05130316	2:3=0.00000000	2:4=0.00000000	3:4=0.00000000
     """
     pos,fsth=parseline(l.rstrip("\n"))
     positive=False
     if(pos%100==1 and pos > 1):
          positive=True
          positivesize+=1
     else:
          negativesize+=1
     
     for key in comp:
          val=fsth[key]
  
          if positive:
               tp.append(val)
          else:
               fp.append(val)

fp=sorted(fp,key=lambda i:-i)
tp=sorted(tp,key=lambda i:-i)

tptp=[]
tpfp=[]
possteps=parse_steps(args.fractions,positivesize)
negsteps=parse_steps(args.fractions,negativesize)
for i,val in enumerate(tp):
     if (i+1) in possteps:
          tptp.insert(0,val)

for i,val in enumerate(fp):
     if (i+1) in negsteps:
          tpfp.insert(0,val)

tptp.insert(0,len(tp))
tpfp.insert(0,len(fp))

ps="\t".join([str(t) for t in reversed(sorted(possteps))])
ns="\t".join([str(t) for t in reversed(sorted(negsteps))])
part2="\t".join([str(t) for t in tptp])
part3="\t".join([str(t) for t in tpfp]) 
print ps+"\t"+part2 +"\t"+ns+"\t"+ part3     
