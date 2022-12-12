import sys
import re


id="b:trala{d:{c:tralal}}"

m=re.search(r"^([^:]*):(.*)$",id)  ### RECURSION SAVE; SPLIT only at first occurence
print m.group(1)
print m.group(2)

"""

l=list("abcdefghij")
del l[9]
print l
a,b=re.split(r"=",s)
print(a,b)
if b.startswith("$"):
    print b[1:]
s="$bla"

m=re.search(r"([+-])",s)
if m is None:
    print ("none")

strand=m.group(1)
print(strand)


tsd=4




tosort=[(1,"bla"),(10,"bla"),(33,"tra"),(5,"hello")]
s=sorted(tosort,key=lambda i:-i[0])
print(s)



print (insert(ref,toins,5,0))
print (insert(ref,toins,5,3))
print (insert(ref,toins,5,4))



def insert(ref,what,pos,tsd):
    left=ref[:pos]
    tlp=pos-int((tsd+1)/2)
    trp=tlp+tsd
    tsdseq=ref[tlp:trp]
    right=ref[pos:]
    return left+tsdseq+what+tsdseq+right

#!/usr/bin/env python
from scipy.stats import chisquare,chi2_contingency
import numpy as np

#print(chisquare([13, 12], f_exp=[12,12]))

#print(chisquare([16, 18, 16, 14, 12, 12], f_exp=np.array([16, 16, 16, 16, 16, 8])))
#print(chisquare([135,2],f_exp=np.array([6,1])))
print(chi2_contingency(np.array([[135, 2], [6,1]])))
#print(np.array([[135, 2], [6,1]]))

print(chi2_contingency(np.array([[135, 2], [5,1]])))

# results (0.040000000000000001, 0.84148058112179402)
# R chisq.test(c(13,12)) =>  X-squared = 0.04, df = 1, p-value = 0.8415
# chisq.test(c(13,12),p=c(12,12),rescale.p=TRUE) => X-squared = 0.04, df = 1, p-value = 0.8415
"""