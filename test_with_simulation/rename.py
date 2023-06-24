import glob
from subprocess import Popen, PIPE, STDOUT, DEVNULL

# rename tp_f_oid.bam --> fp_oid.bam, with oids in nohup.out
oids = set()
for l in open("nohup.out", "r"):
    if l.startswith('all'):
        oid      = l.split()[4]
        # cmd      = ['mv tp_f_{}.bam fp_{}.bam && mv tp_f_{}.bam.bai fp_{}.bam.bai'.format(oid, oid, oid, oid)]
        # idx_proc = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
        # exitcode = idx_proc.wait()

        # if exitcode != 0:
        #     raise Exception("Error: rename for {} failed".format(oid))
        oids.add(oid)


# move record in TP.bed into FP.bed, with record's oid in nohup.out
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