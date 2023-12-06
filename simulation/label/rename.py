import glob
from subprocess import Popen, PIPE, STDOUT, DEVNULL

### rename tp_f_oid.bam --> fp_oid.bam, with oids in log.txt ###
oids = set()
try:
    for l in open("log.txt", "r"):
        if l.startswith('all'):
            oid = l.split()[4]
            oids.add(oid)
except:
    print("no log.txt")

### move record in TP.bed into FP.bed, with record's oid in log.txt ###
tpnew = open("TP_new.bed", "w")
fp    = open("FP.bed", "a")

for l in open("TP.bed", "r"):
    l   = l.strip().split()
    oid = l[8]
    if oid in oids:
        # write into FP.bed
        fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(l[0],l[6],l[7],
                                                             l[8],l[4],l[5],
                                                             l[1],l[2],l[3]))
    else:
        # write into TP_new.bed
        tpnew.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(l[0],l[1],l[2],
                                                                l[3],l[4],l[5],
                                                                l[6],l[7],l[8]))

tpnew.close()
fp.close()


### split TP & FP cluster records ###
oids = set()
for l in open("TP_new.bed", "r"):
    oid = l.split()[8]
    oids.add(oid)

tp_out = open("TP_clt.txt", "w")
fp_out = open("FP_clt.txt", "w")
for l in open("tmp_all_clt.txt", "r"):
    oid = l.split()[3]
    if oid in oids:
        tp_out.write(l)
    else:
        fp_out.write(l)

tp_out.close()
fp_out.close()
