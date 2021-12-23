#!/bin/bash

chromosomes=$(samtools view -H "${1}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

for chr in $chromosomes; do
    samtools view -@ 48 "${1}" "${chr}" -o "tmp_${chr}."sam
done