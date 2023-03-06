#!/bin/bash

# extract header
samtools view -H all_candidates_alignments.bam > header


# merge fp support alignment
for bam in fp_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o fp_merge.bam -
samtools sort -o tmp.bam fp_merge.bam && mv tmp.bam fp_merge.bam && samtools index fp_merge.bam && rm tmp.sam


# merge true tp support alignment
for bam in tp_t_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o tp_t_merge.bam -
samtools sort -o tmp.bam tp_t_merge.bam && mv tmp.bam tp_t_merge.bam && samtools index tp_t_merge.bam && rm tmp.sam


# merge false tp support alignment
for bam in tp_f_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o tp_f_merge.bam -
samtools sort -o tmp.bam tp_f_merge.bam && mv tmp.bam tp_f_merge.bam && samtools index tp_f_merge.bam && rm tmp.sam