## 1. Generate raw reads

### 1.1 First generate population genome definition files on local machine by running `simulation_protocol.sh`

```shell
# 1. For Testing Germline Insertion
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_germ_ccs
simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/benchmark/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 2500 --avg-somatic-count 0 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2

# 2. For Testing Somatic Insertion
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_soma_ccs
simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/benchmark/GRCh38.transposon.fa -N 1000 -R 0 --sub-N 50 --germline-count 0 --avg-somatic-count 10 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2
```

### 1.2 Then build population genome and generate raw reads by submitting `simulation_with_slurm_GRCh38.sh`

```shell
sbatch simulation_with_slurm_GRCh38.sh
```

### 1.3 Then merge raw reads from each chromosome by submitting `merge_with_slurm_GRCh38.sh`

```shell
sbatch merge_with_slurm_GRCh38.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_germ_ccs:
# 1. TGS.fasta
# 2. NGS_1.fastq
# 3. NGS_2.fastq
# 4. chr1/chr1.ins.summary
# 5. chr1/chr1.0.pgd
```

3 datasets was generated for testing germline insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 100                    | 2500                     | 0                       | 0                         | 10% (0.1~1); 90%(0.5 or 1)  | 0                           | 0.7 (only L1)                  | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

3 datasets was generated for testing somatic insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 1000                   | 0                        | 10000                   | 10                        | 0                           | 0.001                       | 0.7 (only L1)                  | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

## 2. Labeling

This step aims to find all insertion candidates and label them based on ground truth.

### 2.1 First map raw reads to reference genome with minimap2 on blades

```shell
sbtach minimap2_with_slurm_GRCh38.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_germ_ccs/map-hifi:
# 1. TGS.bam
# 2. TGS.bam.bai
```

### 2.2 Find all candidates on local machine

```shell
bash batch_labeling_GRCh38.sh
```

### 2.3 Label all candidates on blades

```shell
# For germline(high frequency) insertion datasets
sbatch label_with_slurm_germline_GRCh38.sh

# For somatic(low frequency) insertion datasets
sbatch label_with_slurm_somatic_GRCh38.sh

# Main results under /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_germ_ccs/map-hifi/result_no_secondary:
# 1. TP_clt.txt
# 2. FP_clt.txt
```
