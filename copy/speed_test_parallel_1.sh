#!/bin/bash

echo -e "speed test with filterring for '${1}': "
start=$(date +%s)
# Get BAM file header
samtools view -H ${1} > tmp.header
time samtools view -@ 48 ${1} > tmp.sam

# split input file
time parallel --pipepart -a tmp.sam --block 10M -j48 -q awk '{split($6, a, /[ISH]/); max=0; for(i in a){match(a[i], /([0-9]+$)/, b); if(b[1]>max){max=b[1]}}; if(max>=200){print $0}}' | cat tmp.header - | samtools view -@ 48 -bhS - | samtools sort -@ 48 -o ${1%.bam}_filtered.bam - && samtools index -@ 48 ${1%.bam}_filtered.bam

python speed_test.py ${1%.bam}_filtered.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds.
# time samtools view -@ 32 test_10t.bam | awk '{split($6, a, /[ISH]/); max=0; for(i in a){match(a[i], /([0-9]+$)/, b); if(b[1]>max){max=b[1]}}; if(max>200){print $0}}' - | cat tmp.header - | samtools view -@ 32 -bhS - | samtools sort -@ 32 -o test_10t_filtered.bam - && samtools index -@ 32 test_10t_filtered.bam


################################################################################
echo -e "speed test without filterring for '${1}': "
start=$(date +%s)

python speed_test.py ${1}

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds.

rm *tmp*
rm *filter*