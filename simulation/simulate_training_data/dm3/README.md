## 1. Generate raw reads

### 1.1 First generate population genome definition files on local machine by running `simulation_protocol.sh`

```shell
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/v4/line_28_0_ccs
/data/tusers/zhongrenhu/for_SMS/bin/simulation/simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 50 --min-distance 9000 --depth 50 --species fly --protocol ccs --truncProb 0.1 --nestProb 0.1 --mode 1
```

### 1.2 Then build population genome and generate raw reads by submitting `simulation_with_slurm_dm3.sh`

```shell
sbatch simulation_with_slurm_dm3.sh
```

### 1.3 Then merge raw reads from each chromosome by submitting `merge_with_slurm_dm3.sh`

```shell
sbatch merge_with_slurm_dm3.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/v4/line_28_0_ccs:
# 1. TGS.fasta
# 2. chr2L/chr2L.ins.summary
# 3. chr2L/chr2L.0.pgd
```

300 datasets was generated for Drosophila using different `depth(10-50X)` and `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth  | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ------ | ------------ | -------------------------------- | --------------------------- |
| 100                    | 1000                     | 5000                    | 50                        | 0.1~1                       | 0.01                        | 0.1                            | 0.1                          | 0                         | 6~8        | 10~50X | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

## 2. Labeling

This step aims to find all insertion candidates and label them based on ground truth.

### 2.1 First map raw reads to reference genome with minimap2 on blades

```shell
sbtach minimap2_with_slurm_dm3.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/v4/line_28_0_ccs/map-hifi:
# 1. TGS.bam
# 2. TGS.bam.bai
```

### 2.2 First find all candidates on local machine

```shell
bash batch_labeling_dm3.sh
```

### 2.3 Label all candidates on blades

```shell
sbatch label_with_slurm_dm3.sh
```

### 2.4 Classify germline and somatic insertion

This step aims to distinguish which positive candidates overlap with germline insertions and which ones overlap with somatic insertions.

```shell
bash classify_germline_and_smoatic.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/v4/line_28_0_ccs/map-hifi/result_no_secondary:
# 1. TP_clt_G.txt
# 2. TP_clt_S.txt
# 3. TP_clt.txt
# 4. FP_clt.txt
```
