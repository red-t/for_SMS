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
parser.add_argument("--tomask", type=int, required=False, dest="tomask", default=5, help="The end to mask")
args = parser.parse_args()
tm =args.tomask

for line in args.sam:
     """
     0  1        2      3       4         5     6        7       8
     1	99	2R	1	0	100M	=	201	300
     1	99	2R	1	55	100M	=	175	274	CTCAAGATACCTTCTACAGATTATTTAAAGCTAGTGCACAACAACAATAAATTGACTAAGTTATGTCATTTTAAGCGGTCAAAATGGGTGATTTTTCGAT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	AS:i:0	UQ:i:0	NM:i:0	MD:Z:100	PQ:i:35	SM:i:1	AM:i:1	PG:Z:novoalign

     """
     a=line.rstrip("\n").split("\t")
     read=a[9]
     l=list(read)
     l[0:tm]="N"*tm
     l[-tm:]="N"*tm
     read="".join(l)
     a[9]=read
     print "\t".join(a)
     
     
     
     
     
     
