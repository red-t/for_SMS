#!/bin/bash

frac1=0.000572311   # 10000
frac2=0.00572311    # 100000
frac3=0.0572311     # 1000000
samtools view -H  ../HG002_hs37d5_ONT-UL_GIAB_20200204.bam > tmp.header
samtools view -@ 24 -s $frac1 ../HG002_hs37d5_ONT-UL_GIAB_20200204.bam | cat tmp.header - | samtools view -@ 24 -bhS > test_10t.bam
samtools view -@ 24 -s $frac2 ../HG002_hs37d5_ONT-UL_GIAB_20200204.bam | cat tmp.header - | samtools view -@ 24 -bhS > test_100t.bam
samtools view -@ 24 -s $frac3 ../HG002_hs37d5_ONT-UL_GIAB_20200204.bam | cat tmp.header - | samtools view -@ 24 -bhS > test_1000t.bam





split -a 2 -b 100m test_100t_split.



awk --profile '{cig=$6; while (match(cig, /([0-9]+)[ISH]/, a)){if(a[1]>=2){print $0; break}; cig=substr(cig, RSTART+RLENGTH)}}'