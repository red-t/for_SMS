## 1. For Drosophila

First generate population genome definition files on local machine by running `simulation_protocol.sh`.

```shell
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/v4/line_28_0_ccs

simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 50 --min-distance 9000 --depth 50 --species fly --protocol ccs

# This script will make a directory for each chromosome, like chr2L, which includes:
chr2L.ins.summary				# summary information of each simulated insertion
chr2L.tmp.chasis.fasta			# sequence of chr2L
chr2L.tmp.chasis.fasta.fai
chr2L.0.pgd						# population genome definition file of the first sub-population genome
chr2L.1.pgd
...
```

Then build population genome and generate raw reads by submitting `simulation_with_slurm_dm3.sh`.

```shell
sbatch simulation_with_slurm_dm3.sh

# Each job will generate raw reads from one chromosome, like:
chr2L/TGS.fasta
```

Then merge raw reads from each chromosome by submitting `merge_with_slurm_dm3.sh`.

```shell
sbatch merge_with_slurm_dm3.sh

# Each job will merge all chromosomes' raw reads into a single FASTA file, like:
line_28_0_ccs/TGS.fasta
```

300 datasets was generated for Drosophila using different `depth(10-50X)` and `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth  | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ------ | ------------ | -------------------------------- | --------------------------- |
| 100                    | 1000                     | 5000                    | 50                        | 0.1~1                       | 0.01                        | 0.1                            | 0.1                          | 0                         | 6~8        | 10~50X | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

## 2. For Human

First generate population genome definition files on local machine by running `simulation_protocol.sh`.

```shell
# Current working directory is /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/v4/HG02716_0_ccs

simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 3000 --avg-somatic-count 500 --min-distance 9000 --depth 50 --species human --protocol ccs
```

Then build population genome and generate raw reads by submitting `simulation_with_slurm_GRCh38.sh`.

```shell
sbatch simulation_with_slurm_GRCh38.sh
```

Then merge raw reads from each chromosome by submitting `merge_with_slurm_GRCh38.sh`.

```shell
sbatch merge_with_slurm_GRCh38.sh
```

150 datasets was generated for Human using different `depth(10-50X)` and `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth  | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ------ | ------------ | -------------------------------- | --------------------------- |
| 100                    | 3000                     | 50000                   | 500                       | 10% (0.1~1); 90%(0.5 or 1)  | 0.01                        | 0.1                            | 0.1                          | 0                         | 6~8        | 10~50X | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |
