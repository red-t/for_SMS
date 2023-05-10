# 1. 筛选 “位于边界处的 TP & FP cluster”
# 2. 筛选 “repeats 内部的 TP & FP cluster”
# 3. 筛选 “正常区域的 FP cluster”

tp1 = open("to_check/TP_at_boundary.bed", "w")
tp2 = open("to_check/TP_inside_repeats.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    ctype = l[6]
    if ctype in ['1', '2', '3', '5']:
        tp1.write("\t".join(l) + "\n")
    if ctype in ['4', '5']:
        tp2.write("\t".join(l) + "\n")

tp1.close()
tp2.close()

fp1 = open("to_check/FP_at_boundary.bed", "w")
fp2 = open("to_check/FP_inside_repeats.bed", "w")
fp3 = open("to_check/FP_at_normal.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    ctype = l[6]
    if ctype == '0':
        fp3.write("\t".join(l) + "\n")
    if ctype in ['4', '5']:
        fp2.write("\t".join(l) + "\n")
    if ctype in ['1', '2', '3', '5']:
        fp1.write("\t".join(l) + "\n")

fp1.close()
fp2.close()
fp3.close()

# ------------------------------------------------------------------------
# 4. 筛选 “包含 mid-insert segments 的 TP & FP clusters”
tp1 = open("to_check/TP_contain_mid.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    nmid = int(l[16])
    if nmid > 0:
        tp1.write("\t".join(l) + "\n")

tp1.close()

fp1 = open("to_check/FP_contain_mid.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    nmid = int(l[16])
    if nmid > 0:
        fp1.write("\t".join(l) + "\n")

fp1.close()

# ------------------------------------------------------------------------
# 5. 筛选 “类型单一的 TP clusters“
# 6. 筛选 “多种类型的 FP clusters“
tp1 = open("to_check/TP_single_type.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    ntype = int(l[12])
    if ntype == 1:
        tp1.write("\t".join(l) + "\n")

tp1.close()

fp1 = open("to_check/FP_multi_type.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    ntype = int(l[12])
    if ntype > 1:
        fp1.write("\t".join(l) + "\n")

fp1.close()

# ------------------------------------------------------------------------
# 7. 筛选 “ratio(<0.51, f1 peak) 较小的 TP“
# 8. 筛选 “ratio(>=0.51, f1 peak) 较大的 FP“
tp1 = open("to_check/TP_low_ratio.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    ratio = float(l[14])
    if ratio < 0.51:
        tp1.write("\t".join(l) + "\n")

tp1.close()

fp1 = open("to_check/FP_high_ratio.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    ratio = float(l[14])
    if ratio >= 0.51:
        fp1.write("\t".join(l) + "\n")

fp1.close()

# ------------------------------------------------------------------------
# 9. 筛选 “entropy 较小的 (<0.95, f1 peak) TP”
# 10. 筛选 “entropy 较大的 (>=0.95, f1 peak) FP”
tp1 = open("to_check/TP_low_entropy.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    entropy = float(l[13])
    if entropy < 0.95:
        tp1.write("\t".join(l) + "\n")

tp1.close()

fp1 = open("to_check/FP_high_entropy.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    entropy = float(l[13])
    if entropy >= 0.95:
        fp1.write("\t".join(l) + "\n")

fp1.close()

# ------------------------------------------------------------------------
# 11. 筛选出 “avg_mapq 较低 (<57)” 的 TP cluster
# 12. 筛选出 “avg_mapq 较高 (>=57)” 的 FP cluster
tp1 = open("to_check/TP_low_mapq.bed", "w")
for l in open("TP_stats.bed", "r"):
    l = l.split()
    avg_mapq = float(l[20])
    if avg_mapq < 57:
        tp1.write("\t".join(l) + "\n")

tp1.close()

fp1 = open("to_check/FP_high_mapq.bed", "w")
for l in open("FP_stats.bed", "r"):
    l = l.split()
    avg_mapq = float(l[20])
    if avg_mapq >= 57:
        fp1.write("\t".join(l) + "\n")

fp1.close()