from subprocess import Popen, PIPE, STDOUT, DEVNULL
import ailist as aiL
import math

# --------------------------------------------------------------------------------------

def parse_stype(stype):
    ret = ''
    if stype & 0x1:
        ret = 'L'
    if stype & 0x2:
        ret = 'M'
    if stype & 0x4:
        ret = 'R'
    return ret

def loc_type(n1, l1, n2, l2):
    if n1==0 and n2==0:
        # normal
        t = 0
    elif (n1>0 and l1<=50) and n2==0:
        # repeat edge
        t = 1
    elif (n2>0 and l2<=50) and n1==0:
        # gap edge
        t = 2
    elif (n1>0 and l1<=50) and (n2>0 and l2<=50):
        # repeat & gap edge
        t = 3
    elif (n1>0 and l1>50) and n2==0:
        # inside repeat
        t = 4
    elif (n1>0 and l1>50) and (n2>0 and l2<=50):
        # gap edge & inside repeat
        t = 5
    return str(t)

def judge_type(c):
    segn = int(c[4])
    ntype = 0
    entropy = 0

    # compute entropy & ntype
    for i in range(15,18):
        if c[i] != 0:
            entropy -= (c[i]/segn) * math.log2(c[i]/segn)                 # compute entropy
            ntype   += 1                                                  # compute ntype

    # compute balance_ratio
    balance_ratio = (min(c[15], c[17])+0.01)/(max(c[15], c[17])+0.01)
    
    # return ntype, entropy, balance_ratio
    return str(ntype), str(entropy), str(balance_ratio)

# --------------------------------------------------------------------------------------

# 1. fill TP with cid of true positive clusters
TP = set()
for l in open("TP.bed", "r"):
    l = l.strip().split()
    TP.add(l[-1])

# 2. construct AIList with repeat/gap
repeat = aiL.ailist()
repeat.addAll('repeat_200.bed')
repeat.construct()
gap = aiL.ailist()
gap.addAll('dm3_gap.bed')
gap.construct()

# 3. to calculate the distance between segment breakpoint and repeat boundary
tout = open("TP_seg.txt", "w")
fout = open("FP_seg.txt", "w")
for l in open("all_candidates_segments.txt", "r"):
    # parse segment line
    l     = l.split()
    chr   = l[0]
    flag  = l[1]
    mapq  = l[2]
    qlen  = int(l[4])-int(l[3])
    rpos  = int(l[5])
    stype = int(l[6])
    cid   = l[7]
    rst   = int(l[11])
    red   = int(l[12])
    qname = l[-13]; NM = l[-8]; AS = l[-7]; ms = l[-6]
    cm    = l[-5];  s1 = l[-4]; rl = l[-2]; de = l[-1]

    # compute distance of rpos to the nearest repeat/gap
    n1, l1 = repeat.dist_p(chr, rpos, 50)
    n2, l2 = gap.dist_p(chr, rpos, 50)

    # compute distance of rst/red to the nearest repeat/gap
    if stype & 0x1:
        # left-clip segment
        n3 = n1; l3 = l1
        n5 = n2; l5 = l2
        n4, l4 = repeat.dist_p(chr, red, 50)
        n6, l6 = gap.dist_p(chr, red, 50)
    elif stype & 0x4:
        # right-clip segment
        n4 = n1; l4 = l1
        n6 = n2; l6 = l2
        n3, l3 = repeat.dist_p(chr, rst, 50)
        n5, l5 = gap.dist_p(chr, rst, 50)
    else:
        # mid-insert segment
        n3, l3 = repeat.dist_p(chr, rst, 50)
        n4, l4 = repeat.dist_p(chr, red, 50)
        n5, l5 = gap.dist_p(chr, rst, 50)
        n6, l6 = gap.dist_p(chr, red, 50)
    
    # transform stype(int) --> stype(str)
    stype = parse_stype(stype)

    # write out
    if cid in TP:
        tout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'
                   '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    chr, flag, mapq, rpos, stype, cid, n1, l1, n2, l2, qname, rst, 
                    red, NM, AS, ms, cm, s1, rl, de, n3, l3, n4, l4, n5, l5, n6, l6, qlen))
    else:
        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'
                   '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    chr, flag, mapq, rpos, stype, cid, n1, l1, n2, l2, qname, rst, 
                    red, NM, AS, ms, cm, s1, rl, de, n3, l3, n4, l4, n5, l5, n6, l6, qlen))

tout.close()
fout.close()

# --------------------------------------------------------------------------------------

# 1. initialize highfreq
highfreq = set()
for l in open("high_freq_ins.bed", "r"):
    l = l.split()
    id = l[0]+l[3]
    highfreq.add(id)

# 2. initialize TP
TP = dict()
for l in open("TP.bed", "r"):
    l   = l.strip().split()
    chr = l[0]
    st  = l[1]
    ed  = l[2]
    id  = l[3]
    nseg   = l[4]
    strand = l[5]
    ost = l[6]
    oed = l[7]
    cid = l[8]

    TP[cid] = [0]*21
    TP[cid][0] = chr
    TP[cid][1] = ost
    TP[cid][2] = oed
    TP[cid][3] = cid
    TP[cid][4] = nseg
    TP[cid][5] = strand
    TP[cid][19] = []
    TP[cid][20] = id

    # judge high/low freq
    t = chr + id
    if t in highfreq:
        TP[cid][11] = 'H'
    else:
        TP[cid][11] = 'L'

# 3. calculate distance for TP clusters
del highfreq
for cid in TP:
    chr = TP[cid][0]
    ost = int(TP[cid][1])
    oed = int(TP[cid][2])

    n1, l1 = repeat.dist_c(chr, ost, oed, 50)
    n2, l2 = gap.dist_c(chr, ost, oed, 50)

    TP[cid][6] = loc_type(n1, l1, n2, l2)
    TP[cid][7] = n1
    TP[cid][8] = l1
    TP[cid][9] = n2
    TP[cid][10] = l2

# 4. read TP segment records
for l in open("TP_seg.txt", "r"):
    l     = l.split()
    mapq  = int(l[2])
    stype = l[4]
    cid   = l[5]
    qlen  = int(l[-1])
    
    # parse segment type
    if stype == 'L':
        TP[cid][15] += 1
    elif stype == 'R':
        TP[cid][17] += 1
    else:
        TP[cid][16] += 1
        TP[cid][19].append(qlen)

    # parse segment mapq
    if mapq<45:
        TP[cid][18] += 1

# 5. writting TP
out = open("TP_stats.bed", "w")
for cid in TP:
    TP[cid]

    # ntype, entropy, balance_ratio
    TP[cid][12], TP[cid][13], TP[cid][14] = judge_type(TP[cid])

    # compute low mapq segment fraction
    TP[cid][18] = str(float(TP[cid][18])/int(TP[cid][4]))

    # compute average length of mid-insert segments
    if(len(TP[cid][19])):
        TP[cid][19] = str(sum(TP[cid][19])/len(TP[cid][19]))
    else:
        TP[cid][19] = str(0)

    # transfer int --> str
    TP[cid][7] = str(TP[cid][7])
    TP[cid][8] = str(TP[cid][8])
    TP[cid][9] = str(TP[cid][9])
    TP[cid][10] = str(TP[cid][10])
    TP[cid][15] = str(TP[cid][15])
    TP[cid][16] = str(TP[cid][16])
    TP[cid][17] = str(TP[cid][17])

    # writting
    out.write("\t".join(TP[cid])+"\n")

out.close()

# --------------------------------------------------------------------------------------

# 1. initialize FP
FP = dict()
for l in open("FP.bed", "r"):
    l   = l.strip().split()
    chr = l[0]
    st  = l[1]
    ed  = l[2]
    cid  = l[3]
    nseg   = l[4]
    strand = l[5]

    FP[cid] = [0]*20
    FP[cid][0] = chr
    FP[cid][1] = st
    FP[cid][2] = ed
    FP[cid][3] = cid
    FP[cid][4] = nseg
    FP[cid][5] = strand
    FP[cid][19] = []

    # judge high/low freq
    FP[cid][11] = 'N'

# 3. calculate distance for FP clusters
for cid in FP:
    chr = FP[cid][0]
    st = int(FP[cid][1])
    ed = int(FP[cid][2])

    n1, l1 = repeat.dist_c(chr, st, ed, 50)
    n2, l2 = gap.dist_c(chr, st, ed, 50)

    FP[cid][6] = loc_type(n1, l1, n2, l2)
    FP[cid][7] = n1
    FP[cid][8] = l1
    FP[cid][9] = n2
    FP[cid][10] = l2

# 4. read FP segment records
for l in open("FP_seg.txt", "r"):
    l     = l.split()
    mapq  = int(l[2])
    stype = l[4]
    cid   = l[5]
    qlen  = int(l[-1])
    
    # parse segment type
    if stype == 'L':
        FP[cid][15] += 1
    elif stype == 'R':
        FP[cid][17] += 1
    else:
        FP[cid][16] += 1
        FP[cid][19].append(qlen)

    # parse segment mapq
    if mapq<45:
        FP[cid][18] += 1

# 5. writting FP
out = open("FP_stats.bed", "w")
for cid in FP:
    FP[cid]

    # ntype, entropy, balance_ratio
    FP[cid][12], FP[cid][13], FP[cid][14] = judge_type(FP[cid])

    # compute low mapq segment fraction
    FP[cid][18] = str(float(FP[cid][18])/int(FP[cid][4]))

    # compute average length of mid-insert segments
    if(len(FP[cid][19])):
        FP[cid][19] = str(sum(FP[cid][19])/len(FP[cid][19]))
    else:
        FP[cid][19] = str(0)

    # transfer int --> str
    FP[cid][7] = str(FP[cid][7])
    FP[cid][8] = str(FP[cid][8])
    FP[cid][9] = str(FP[cid][9])
    FP[cid][10] = str(FP[cid][10])
    FP[cid][15] = str(FP[cid][15])
    FP[cid][16] = str(FP[cid][16])
    FP[cid][17] = str(FP[cid][17])

    # writting
    out.write("\t".join(FP[cid])+"\n")

out.close()

# --------------------------------------------------------------------------------------







# calculate distance for TP clusters
# fpout = open("TP_with_dist.bed", "w")
# for l in open("TP.bed", "r"):
#     l   = l.strip().split()
#     chr = l[0]
#     st  = int(l[6])
#     ed  = int(l[7])
#     cid = l[3]
#     nseg   = l[4]
#     strand = l[5]

#     n1, l1 = repeat.dist_c(chr, st, ed, 50)
#     n2, l2 = gap.dist_c(chr, st, ed, 50)

#     fpout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
#                 chr, st, ed, cid, nseg, strand, n1, l1, n2, l2))

# fpout.close()

# calculate distance for FP clusters
# fpout = open("FP_with_dist.bed", "w")
# for l in open("FP.bed", "r"):
#     l   = l.strip().split()
#     chr = l[0]
#     st  = int(l[1])
#     ed  = int(l[2])
#     cid = l[3]
#     nseg   = l[4]
#     strand = l[5]

#     n1, l1 = repeat.dist_c(chr, st, ed, 50)
#     n2, l2 = gap.dist_c(chr, st, ed, 50)

#     fpout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
#                 chr, st, ed, cid, nseg, strand, n1, l1, n2, l2))

# fpout.close()
# del repeat; del gap

# --------------------------------------------------------------------------------------

# high = set()
# for l in open("high_freq_ins.bed", "r"):
#     l = l.split()
#     id = l[0]+l[3]
#     high.add(id)

# hout = open("htmp", "w")
# lout = open("ltmp", "w")
# for l in open("TP_with_dist.bed"):
#     l = l.strip().split()
#     id = l[0]+l[3]
#     if id in high:
#         hout.write("\t".join(l)+"\n")
#     else:
#         lout.write("\t".join(l)+"\n")

# hout.close()
# lout.close()

# --------------------------------------------------------------------------------------

# tp = dict()
# for l in open('TP_seg.txt', 'r'):
#     l     = l.split()
#     stype = l[4]
#     cid   = l[5]
#     try:
#         tp[cid][0] += 1
#         if stype == "L":
#             tp[cid][1] += 1
#         elif stype == "R":
#             tp[cid][3] += 1
#         else:
#             tp[cid][2] += 1
#     except KeyError:
#         tp[cid] = [1,0,0,0]
#         if stype == "L":
#             tp[cid][1] += 1
#         elif stype == "R":
#             tp[cid][3] += 1
#         else:
#             tp[cid][2] += 1

# out = open("TP_stype_distrb", "w")
# for cid in tp:
#     clt = tp[cid]
#     out.write('{}\t{}\t{}\t{}\t{}\n'.format(cid, clt[0], clt[1], clt[2], clt[3]))

# out.close()

# --------------------------------------------------------------------------------------

# tp = dict()
# for l in open('TP_seg.txt', 'r'):
#     l    = l.split()
#     mapq = int(l[2])
#     cid  = l[5]
#     try:
#         tp[cid][0] += 1
#         if mapq < 30:
#             tp[cid][1] += 1
#         else:
#             tp[cid][2] += 1
#     except KeyError:
#         tp[cid] = [1,0,0]
#         if mapq < 30:
#             tp[cid][1] += 1
#         else:
#             tp[cid][2] += 1

# out = open("TP_mapq_distrb", "w")
# for cid in tp:
#     clt = tp[cid]
#     out.write('{}\t{}\t{}\t{}\n'.format(cid, clt[0], clt[1], clt[2]))

# out.close()

# --------------------------------------------------------------------------------------

# a = []
# a.append('($21==0 && $25==0)')              # st 处于正常区域
# a.append('($21>0 && $22<=50 && $25==0)')    # st 仅位于 repeat 边界
# a.append('($21==0 && $25>0)')               # st 仅位于 gap 边界
# a.append('($21>0 && $22<=50 && $25>0)')     # st 同时位于 repeat 及 gap 边界
# a.append('($21>0 && $22>50 && $25==0)')     # st 仅位于 repeat 内
# a.append('($21>0 && $22>50 && $25>0)')      # st 同时位于 repeat 内及 gap 边界

# b = []
# b.append('($23==0 && $27==0)')              # ed 处于正常区域
# b.append('($23>0 && $24<=50 && $27==0)')    # ed 仅位于 repeat 边界
# b.append('($23==0 && $27>0)')               # ed 仅位于 gap 边界
# b.append('($23>0 && $24<=50 && $27>0)')     # ed 同时位于 repeat 及 gap 边界
# b.append('($23>0 && $24>50 && $27==0)')     # ed 仅位于 repeat 内
# b.append('($23>0 && $24>50 && $27>0)')      # ed 同时位于 repeat 内及 gap 边界

# for i in range(len(a)):
#     for j in range(len(b)):
#         judg = "{}&&{}".format(a[i], b[j])
#         outf = "tmp/FP_{}-{}".format(i,j)
#         cmd  = ["awk '" + judg + "{print $1,$11,$12,$3,$14,$15,$16,$17,$18,$19,$20}' FP_seg.txt >> " + outf]
#         proc = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
#         exitcode = proc.wait()

#         cmd  = ["sort -u {} > tt && mv tt {}".format(outf, outf)]
#         proc = Popen(cmd, stderr=DEVNULL, shell=True, executable='/bin/bash')
#         exitcode = proc.wait()

# fout = open("TP_alnedge_distrb", "w")
# for i in range(6):
#     for j in range(6):
#         id = '{}-{}'.format(i,j)
#         n  = 0
#         nlow = 0
#         nhigh = 0
#         fn = 'tmp/TP_{}-{}'.format(i,j)
#         for l in open(fn, "r"):
#             if l:
#                 l = l.split()
#                 n += 1
#                 mapq = int(l[3])
#                 if mapq<30:
#                     nlow += 1
#                 else:
#                     nhigh+= 1
        
#         fout.write('{}\t{}\t{}\t{}\n'.format(id, n, nlow, nhigh))

# fout.close()