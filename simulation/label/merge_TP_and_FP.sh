#!/bin/bash

# extract header
samtools view -H $1 > header


# merge fp support alignment
for bam in fp_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o merge_fp.bam -
samtools sort -o tmp.bam merge_fp.bam && mv tmp.bam merge_fp.bam && samtools index merge_fp.bam && rm tmp.sam


# merge true tp support alignment
for bam in tp_t_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o merge_tp_t.bam -
samtools sort -o tmp.bam merge_tp_t.bam && mv tmp.bam merge_tp_t.bam && samtools index merge_tp_t.bam && rm tmp.sam


# merge false tp support alignment
for bam in tp_f_*bam
do
    samtools view $bam >> tmp.sam
done

cat header tmp.sam | samtools view -b -o merge_tp_f.bam -
samtools sort -o tmp.bam merge_tp_f.bam && mv tmp.bam merge_tp_f.bam && samtools index merge_tp_f.bam && rm tmp.sam

if [ -f 'TP_new.bed' ];then
    mv TP_new.bed TP.bed
fi