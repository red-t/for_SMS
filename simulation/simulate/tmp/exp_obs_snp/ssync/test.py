#!/usr/bin/env python
import re

cig="3M31I15S18D3M"

prog = re.compile(r"(\d+)([MISDHN])")
result = re.findall(prog,cig)

bla=0;




