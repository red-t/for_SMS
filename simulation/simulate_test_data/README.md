## 1. For Drosophila

First generate population genome definition files on local machine by running `simulation_protocol.sh`.

```shell
# 1. For Testing Germline Insertion
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/benchmark/line_28_germ_ccs
simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/benchmark/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 500 --avg-somatic-count 0 --min-distance 9000 --depth 50 --species fly --protocol ccs --truncProb 0.1 --nestProb 0 --mode 2

# 2. For Testing Somatic Insertion
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/benchmark/line_28_soma_ccs
simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/benchmark/dm3.transposon_for_simulaTE.fa -N 1000 -R 0 --sub-N 50 --germline-count 0 --avg-somatic-count 10 --min-distance 9000 --depth 50 --species fly --protocol ccs --truncProb 0.1 --nestProb 0 --mode 2
```

Then build population genome and generate raw reads by submitting `simulation_with_slurm_dm3.sh`.

```shell
sbatch simulation_with_slurm_dm3.sh
```

Then merge raw reads from each chromosome by submitting `merge_with_slurm_dm3.sh`.

```shell
sbatch merge_with_slurm_dm3.sh
```

3 datasets was generated for testing germline insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 100                    | 500                      | 0                       | 0                         | 0.1~1                       | 0                           | 0.1                            | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

3 datasets was generated for testing somatic insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 1000                   | 0                        | 10000                   | 10                        | 0                           | 0.001                       | 0.1                            | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

## 2. For Human

First generate population genome definition files on local machine by running `simulation_protocol.sh`.

```shell
# 1. For Testing Germline Insertion
# Current working directory is /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_germ_ccs
simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/benchmark/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 2500 --avg-somatic-count 0 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2

# 2. For Testing Somatic Insertion
# Current working directory is /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/HG02716_soma_ccs
simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/benchmark/GRCh38.transposon.fa -N 1000 -R 0 --sub-N 50 --germline-count 0 --avg-somatic-count 10 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2
```

Then build population genome and generate raw reads by submitting `simulation_with_slurm_GRCh38.sh`.

```shell
sbatch simulation_with_slurm_GRCh38.sh
```

Then merge raw reads from each chromosome by submitting `merge_with_slurm_GRCh38.sh`.

```shell
sbatch merge_with_slurm_GRCh38.sh
```

3 datasets was generated for testing germline insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 100                    | 2500                     | 0                       | 0                         | 10% (0.1~1); 90%(0.5 or 1)  | 0                           | 0.7 (only L1)                  | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

3 datasets was generated for testing somatic insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 1000                   | 0                        | 10000                   | 10                        | 0                           | 0.001                       | 0.7 (only L1)                  | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |
