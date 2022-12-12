import sys
import random
import re

# Transfer degenerate bases into normal bases
# python change_degenerate_bases.py input.fa output.fa

fout = open(sys.argv[2], "w")
for l in open(sys.argv[1], "r"):
    if not l.startswith(">"):
        l = l.upper()
        l = re.sub('R', random.choice(['A', 'G']), l)
        l = re.sub('Y', random.choice(['C', 'T']), l)
        l = re.sub('M', random.choice(['A', 'C']), l)
        l = re.sub('K', random.choice(['G', 'T']), l)
        l = re.sub('S', random.choice(['G', 'C']), l)
        l = re.sub('W', random.choice(['A', 'T']), l)
        l = re.sub('H', random.choice(['A', 'T', 'C']), l)
        l = re.sub('B', random.choice(['G', 'T', 'C']), l)
        l = re.sub('V', random.choice(['G', 'A', 'C']), l)
        l = re.sub('D', random.choice(['G', 'A', 'T']), l)
        fout.write(l)
    else:
        fout.write(l)

fout.close()