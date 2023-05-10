----
## For TP & FP filtering

#### 1. run `Filter_TP_and_FP.py`
```
# 通过 “查找重叠区间” + source-target hg 判定，初步筛选出 TP & FP cluster candidates
# 写出相应 segments 的信息
# 写出相应 alignments 的信息

# the working directory should contain:
# 1. all_candidates_alignments.bam
# 2. all_candidates_alignments.bam.bai
# 3. all_candidates.bed
# 4. AllIns.bed.gz
# 5. AllIns.bed.gz.tbi
# 6. all_candidates_segments.txt

nohup python Filter_TP_and_FP.py &
```

#### 2. run `rename.py`
```
# 部分初步筛选得到的 TP，当中包含的 support segments 全部被判定为 “假的(False) support segments”，这部分 TP 实际上是 FP
# 将这部分 TP 修改为 FP (包括相应的 alignments)

# records in oid in nohup.out have no "true support alignments", should be false postive
# 1. to rename tp_f_oid.bam --> fp_oid.bam, with oids in nohup.out
# 2. move record in TP.bed into FP.bed, with record's oid in nohup.out

nohup python rename.py &

# result should include:
# 1. all_candidates.bed
# 2. FP.bed
# 3. TP.bed
```

#### 3. merge results
```
# 将每个 TP & FP 对应的 alignments (BAM) 进行合并，方便查看

nohup ./merge_TP_and_FP.sh

# results should include:
# 1. fp_merge.bam
# 2. fp_merge.bam.bai
# 3. tp_f_merge.bam
# 4. tp_f_merge.bam.bai
# 5. tp_t_merge.bam
# 6. tp_t_merge.bam.bai
```

----

## For Sup & Uns filtering

#### 1. run `Filter_Sup_and_Uns.py`
```
nohup python Filter_Sup_and_Uns.py &
```

#### 2. merge results
```
nohup ./merge_Sup_and_Uns.sh

# results should contain:
# 1. sup_merge.bam
# 2. sup_merge.bam.bai
# 3. uns_merge.bam
# 4. uns_merge.bam.bai
```

----

## For statistics calculating

#### 1. Compile AIList
```
# 需要计算的 特征值(或统计量) 涉及查找重叠区间，因此直接修改 AIList 拿来用，后面整合的时候也更容易

# the working directory should contain:
# 1. AIList/ail.pyx
# 2. AIList/AIList.c
# 3. AIList/AIList.h
# 4. AIList/khash.h
# 5. setup.py

python setup.py build
```

#### 2. run `compute_seg_stats.py`
```
# 为每个 TP & FP segments 计算特征，并将其写出
# 为每个 TP & FP cluster 计算特征，并将其写出

# the working directory should contain:
# 1. TP.bed
# 2. FP.bed
# 3. repeat_200.bed
# 4. dm3_gap.bed
# 5. all_candidates_segments.txt
# 6. high_freq_ins.bed

nohup python compute_seg_stats.py &

# result should include:
# 1. TP_seg.txt
# 2. TP_stats.bed
# 3. FP_seg.txt
# 4. FP_stats.bed

# TP_stats.bed（22列）:
# chr   start   end cid nseg    strand  ctype   n1  l1  n2  l2  freq_type   ntype   entropy balance_ratio   nL  nM  nR  low_mapq_frac   avg_midsegL avg_mapq    id
# FP_stats.bed（21列）:
# chr   start   end cid nseg    strand  ctype   n1  l1  n2  l2  freq_type   ntype   entropy balance_ratio   nL  nM  nR  low_mapq_frac   avg_midsegL avg_mapq
```