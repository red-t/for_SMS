#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
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
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Input
__________
A SNP call file, e.g.:
True	True	2R	201	T	38	1.0	0.5	19:19:0:0
True	True	2R	401	A	38	0.0	0.5	19:0:0:19
True	True	2R	801	A	42	1.0	0.5	21:21:0:0




Output
----------
todo

Authors
-------
    Robert Kofler
""")

class SNP:
     def __init__(self,strbias,reffreq,coverage,A,T,C,G):
          self.strbias=strbias
          self.reffreq=reffreq
          self.coverage=coverage
          self.A=A
          self.T=T
          self.C=C
          self.G=G
          # allele count list=
          self.acl=[('A',A),('T',T),('C',C),('G',G)]
          

def ms(li):
     if(len(li)<1):
          return 0,'na','na'
     mean=0.0
     for s in li:
          mean+=s
     mean=mean/len(li)
     
     sos=0.0
     for s in li:
          difsquare=math.pow(math.fabs(s-mean),2)
          sos+=difsquare
     sos=sos/len(li)
     sos=math.sqrt(sos)
     return len(li), mean, sos
     
def outlier(snplist,maxdiff):
     toret=[]
     lowfreq=0.5-maxdiff
     upfreq=0.5+maxdiff
     for s in snplist:
          if s.reffreq<lowfreq or s.reffreq>upfreq:
               toret.append(s)
     return toret


def filter(snplist,maxbias,mincoverage,maxcoverage=1000000):
     toprocess=[]
     lowbias=0.5-maxbias
     upbias=0.5+maxbias
     for s in snplist:
          if s.coverage<mincoverage:
               continue
          if s.coverage>maxcoverage:
               continue
          if s.strbias<lowbias or s.strbias>upbias:
               continue
          toprocess.append(s)
     
     return toprocess


parser.add_argument("--snp-calls", type=str, required=True, dest="snpcalls", default=None, help="A SNP call file")
args = parser.parse_args()

totalsimulated=19999
totalcount=0
truepositives=0
falsepositives=0
falsenegatives=0
covcla=collections.defaultdict(lambda:0)

tplist=[]
fnlist=[]
fplist=[]

for l in open(args.snpcalls):
     a=l.rstrip("\n").split("\t")
     """
     0            1      2       3      4        5       6       7         8
     True	True	2R	801	A	42	1.0	0.5	21:21:0:0
     """
     simulated=a[0]=="True"
     polymorphic=a[1]=="True"
     pos=int(a[3])
     coverage=int(a[5])
     if(a[6]=="Na" or a[7]=="Na"):
          continue
     ckey=int(coverage/10)
     covcla[ckey]+=1
     strbias=float(a[6])
     reffreq=float(a[7])
     totalcount+=1

     A,T,C,G=[int(i) for i in a[8].split(":")]
     snp=SNP(strbias,reffreq,coverage,A,T,C,G)
     
     # true positives false positives etc
     if simulated:
          if polymorphic:
               truepositives+=1
               tplist.append(snp)
          else:
               falsenegatives+=1
               fnlist.append(snp)
     elif polymorphic:
          falsepositives+=1
          fplist.append(snp)
     else:
          raise Exception("Invalid SNP status"+simulated+polymorphic)
     
     




tmp=[]
for i in range(0,max(covcla.keys())+1):
     if i in covcla:
          toadd=str(covcla[i])
          if i==20:
               toadd="*"+toadd+"*"
          tmp.append(toadd)
     else:
          tmp.append("0")
covstr=" ".join(tmp)

print "SNPs\tsimulated\t{0}".format(totalsimulated)
print "SNPs\tfound(P)\t{0}".format(totalcount)
print "SNPs\ttrue_found(TP)\t{0}".format(truepositives)
print "SNPs\ttrue_missed(FN)\t{0}".format(totalsimulated-truepositives)
print "SNPs\ttrue_missed_snpcall(FN1)\t{0}".format(falsenegatives)
print "SNPs\ttrue_missed_pileup(FN2)\t{0}".format(totalsimulated-truepositives-falsenegatives)
print "SNPs\tfalse_found(FP)\t{0}".format(falsepositives)
c,mrf, stdrf=ms([s.reffreq for s in tplist])
print("TP\tmean_reffreq\t{0}".format(mrf))
print("TP\tstd_reffreq\t{0}".format(stdrf))
    
c,mrf, stdrf=ms([s.strbias for s in tplist])
print("TP\tmean_strbias\t{0}".format(mrf))
print("TP\tstd_strbias\t{0}".format(stdrf))

# outlier

c=len(outlier(tplist,0.05))
print("TP\toutlier(0.05)\t{0}".format(c))    
c=len(outlier(tplist,0.1))
print("TP\toutlier(0.1)\t{0}".format(c))    
c=len(outlier(tplist,0.2))
print("TP\toutlier(0.2)\t{0}".format(c))  
c=len(outlier(tplist,0.3))
print("TP\toutlier(0.3)\t{0}".format(c))  
c=len(outlier(tplist,0.4))
print("TP\toutlier(0.4)\t{0}".format(c))  





c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.3,0)])
print("TP_sb03\tcount\t{0}".format(c))
print("TP_sb03\tfilter_loss\t{0}".format(truepositives-c))
print("TP_sb03\tmean_reffreq\t{0}".format(mrf))
print("TP_sb03\tstd_reffreq\t{0}".format(stdrf))



c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.2,0)])
print("TP_sb02\tcount\t{0}".format(c))
print("TP_sb02\tfilter_loss\t{0}".format(truepositives-c))
print("TP_sb02\tmean_reffreq\t{0}".format(mrf))
print("TP_sb02\tstd_reffreq\t{0}".format(stdrf))


c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.1,0)])
print("TP_sb01\tcount\t{0}".format(c))
print("TP_sb01\tfilter_loss\t{0}".format(truepositives-c))
print("TP_sb01\tmean_reffreq\t{0}".format(mrf))
print("TP_sb01\tstd_reffreq\t{0}".format(stdrf))

c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.05,0)])
print("TP_sb005\tcount\t{0}".format(c))
print("TP_sb005\tfilter_loss\t{0}".format(truepositives-c))
print("TP_sb005\tmean_reffreq\t{0}".format(mrf))
print("TP_sb005\tstd_reffreq\t{0}".format(stdrf))

c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.6,50)])
print("TP_mc50\tcount\t{0}".format(c))
print("TP_mc50\tfilter_loss\t{0}".format(truepositives-c))
print("TP_mc50\tmean_reffreq\t{0}".format(mrf))
print("TP_mc50\tstd_reffreq\t{0}".format(stdrf))

c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.6,100)])
print("TP_mc100\tcount\t{0}".format(c))
print("TP_mc100\tfilter_loss\t{0}".format(truepositives-c))
print("TP_mc100\tmean_reffreq\t{0}".format(mrf))
print("TP_mc100\tstd_reffreq\t{0}".format(stdrf))


c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.6,150)])
print("TP_mc150\tcount\t{0}".format(c))
print("TP_mc150\tfilter_loss\t{0}".format(truepositives-c))
print("TP_mc150\tmean_reffreq\t{0}".format(mrf))
print("TP_mc150\tstd_reffreq\t{0}".format(stdrf))


c,mrf, stdrf=ms([s.reffreq for s in filter(tplist,0.1,150)])
print("TP_sb01_mc150\tcount\t{0}".format(c))
print("TP_sb01_mc150\tfilter_loss\t{0}".format(truepositives-c))
print("TP_sb01_mc150\tmean_reffreq\t{0}".format(mrf))
print("TP_sb01_mc150\tstd_reffreq\t{0}".format(stdrf))

print "SNPs\tcoverage_distribution\t{0}".format(covstr)



    