#!/bin/bash

# extract header
samtools view -H /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28_gamma/line_28_pacbio.bam > header


# merge ture support alignment
for bam in sup_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o sup_merge.bam -
samtools sort -o tmp.bam sup_merge.bam && mv tmp.bam sup_merge.bam && samtools index sup_merge.bam && rm tmp.sam


# merge false support alignment
for bam in uns_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o uns_merge.bam -
samtools sort -o tmp.bam uns_merge.bam && mv tmp.bam uns_merge.bam && samtools index uns_merge.bam && rm tmp.sam