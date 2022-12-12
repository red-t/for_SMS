#!/usr/bin/env python

c="AAATTTTTAAA"

print c

l=list(c)
l[0:3]="N"*3
l[-3:]="N"*3
print l