----
## For TP & FP filtering

#### 1. run `Filter_TP_and_FP.py`
```
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
# records in oid in nohup.out have no "true support alignments", should be false postive
# 1. to rename tp_f_oid.bam --> fp_oid.bam, with oids in nohup.out
# 2. move record in TP.bed into FP.bed, with record's oid in nohup.out

nohup python rename.py &

# result should contain:
# all_candidates.bed
# FP.bed
# TP.bed
# TP_new.bed
```

#### 3. merge results
```
nohup ./merge_TP_and_FP.sh

# results should contain:
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