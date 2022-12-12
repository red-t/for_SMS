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
from scipy.stats import fisher_exact,chi2_contingency,chisquare
import numpy as np

maxpvalue=100



def get_strandbiascomputer(sbalgorithm):
     if sbalgorithm == "maxfraction":
          return get_maxfraction
     elif sbalgorithm == "alleleimbalance":
          return get_alleleimbalance
     elif sbalgorithm == "maxsignificance":
          return get_maxsignificance
     else:
          raise Exception("unknown stranbias algorith "+sbalgorithm)
          

def get_measurecomputer(mealgo):
     if mealgo=="fst":
          return get_fst
     elif mealgo=="fet":
          return get_fet
     else:
          raise Exception("unknown measure calculator "+mealgo)



#######################  STRAND BIAS ##################################

def get_maxsignificance(issnp,pcs):
     maxp=0.0
     for s in pcs.samples:
          fwd=s.get_forward()
          rev=s.get_reverse()
          oddsration,pvalue=chisquare([fwd,rev])
          if not pvalue.__nonzero__():
               return maxpvalue
          pval=(-1)*math.log10(pvalue)
          if pval>maxp:
               maxp=pval
     return maxp

def get_alleleimbalance(issnp,pcs):
     """
       | a1    a2
     + | maj+  min+
     - | maj-  min-
     chisquare([maj+,min+],f_exp=[maj-,min-])
     """
     if not issnp:
          return 0.0
     maj,min=pcs.get_maj_min()
     toret=0.0
     for s in pcs.samples:
          mafwd=s.get_XXX_upper(maj)+1
          mifwd=s.get_XXX_upper(min)+1
          marev=s.get_XXX_lower(maj)+1
          mirev=s.get_XXX_lower(min)+1
          chi,pvalue,tmp,temp = chi2_contingency(np.array([[mafwd,mifwd], [marev,mirev]]))

          
          
          fpval=float(pvalue)
          pval=(-1) * math.log10(pvalue)
          if pval> toret:
               toret=pval
     return toret

def get_maxfraction(issnp,pcs):
     return pcs.get_maxstrandbias()
     


###################### MEASURE #########################################

def get_fet(issnp,pcs):     
     if not issnp:
          return 0.0
     
     if len(pcs.samples)!=2:
          raise Exception("Script only works for two samples")
     # major minor allele
     maj,min=pcs.get_maj_min()
     majc,minc=pcs.get_XXX(maj),pcs.get_XXX(min)
     
     oddsratio, pvalue = fisher_exact([[majc[0], minc[0]], [majc[1], minc[1]]])
     if not pvalue.__nonzero__():
          return maxpvalue
     pval=(-1.0)*math.log10(pvalue)
     
     return pval
     #fst=__get_fst(majc[0],minc[0],majc[1],minc[1])

def get_fst(issnp,pcs):
     if not issnp:
          return 0.0
     
     if len(pcs.samples)!=2:
          raise Exception("Script only works for two samples")
     # major minor allele
     maj,min=pcs.get_maj_min()
     majc,minc=pcs.get_XXX(maj),pcs.get_XXX(min)
     fst=__get_fst(majc[0],minc[0],majc[1],minc[1])
     return fst
     
     
def __get_het(p):
     q=1.0-p
     h=1-p**2-q**2
     return h


def __get_fst(ma1,mi1,ma2,mi2):
     if (ma1+mi1)<0.0000001:
          return 0.0
     if(ma2+mi2)<0.0000001:
          return 0.0
     p1=float(ma1)/float(ma1+mi1)
     p2=float(ma2)/float(ma2+mi2)
     
     h1=__get_het(p1)
     h2=__get_het(p2)
     hbet=(h1+h2)/2.0
     
     p12=(p1+p2)/2.0
     ht=__get_het(p12)
     
     fst=(ht-hbet)/ht
     
     return fst
     
     



parser = argparse.ArgumentParser(description="""           
Possible strandbiases
a.) maxfraction; the maximum strandbias in any of the tested populations; works for SNPs and not SNPs; not considering allele frequency differences

b.) alleleimbalance; only works for SNPs; will report the most significant allelic imbalance in any population

c.) maxsignificance; the maximum significant strandbias in any popopulation;

allelicimbalance:
       | a1    a2
     + | maj+  min+
     - | maj-  min-
     chisquare([maj+,min+],f_exp=[maj-,min-])
highest chisquare for any of the tested populations; 0.0 for non-SNPs

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
parser.add_argument("--min-cov", type=int, required=True, dest="mincov", default=None, help="the minimum coverage in all populations; all other sites are ignored")
parser.add_argument("--sb", type=str, required=True, dest="sb", default="maxfraction", help="which strandbias should be computed [maxfraction, alleleimbalance, maxsignificance]")
parser.add_argument("--measure", type=str, required=True, dest="measure", default="fet", help="which measure to compute [fst, fet]")
parser.add_argument("--min-count",type=int, required=False, dest="mincount", default=2, help="the minimum count; for SNP calling")
args = parser.parse_args()

strandbiascomputer=get_strandbiascomputer(args.sb)
measurecomputer=get_measurecomputer(args.measure)
w=MultimeasureWriter("/dev/stdout")

for s in StrandSyncReader(args.ssync):
     if not s.is_coveragesufficient(args.mincov):
          continue
     chr=s.chr
     pos=s.pos
     issnp=s.isPolymorphic(args.mincount)

     sb=strandbiascomputer(issnp,s)
     measure=measurecomputer(issnp,s)
     m=Multimeasure(chr,pos,issnp,False,sb,measure)
     w.write(m)
w.close()
     


