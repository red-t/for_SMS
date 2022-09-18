#!/bin/bash
start=$(date +%s)

DATA_PATH="/data/tusers/boxu/lrft/rawdata/simulation"
OUT_PATH_dir="/data/tusers/boxu/lrft/result/simulation"
lrft_BIN_PATH="/data/tusers/boxu/lrft/scripts"
ANNO_PATH="/data/tusers/boxu/annotation"
REPEATMASKER_FILE="/data/tusers/boxu/annotation/dm3/dm3.rmsk.bed"
GENOME="dm6_clean"
SAMPLE="line_28"
CPU="10"
meth="lrft"
# python /data/tusers.ds/boxu/lrft2/scripts/setup.py build && cp /data/tusers.ds/boxu/lrft2/scripts/build/lib.linux-x86_64-3.7/read_alignment.cpython-37m-x86_64-linux-gnu.so /data/tusers.ds/boxu/lrft2/bin/
# chr2L:7,501,451-7,502,704

for dep in 10 # 1 2 3 4 5 10 20 30 40 50
do
    # [ -f /data/tusers.ds/boxu/lrft2/bin/time_test_result.txt ] && rm /data/tusers.ds/boxu/lrft2/bin/time_test_result.txt
    
    # break
    start=$(date +%s)
    OUT_PATH=${OUT_PATH_dir}/${dep}
    [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
    # [ ! -d ${OUT_PATH}/tldr ] && mkdir ${OUT_PATH}/tldr
    # [ ! -d ${OUT_PATH}/lrft ] && mkdir ${OUT_PATH}/lrft
    # [ ! -d ${OUT_PATH}/TEMP2 ] && mkdir ${OUT_PATH}/TEMP2
    [ ! -d ${OUT_PATH}/SMS ] && mkdir ${OUT_PATH}/SMS
    [ ! -d ${OUT_PATH}/SMS/temp ] && mkdir ${OUT_PATH}/SMS/temp 
    if [ ! -f ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.sorted.bam ];then
        samtools view -bhSq 0 ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.bam > ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.bam
        samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.sorted.bam ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.bam
        samtools index -@ 10 ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.sorted.bam
    fi

    echo -e "chr2L\t3761809\t3761818" > test.bed
    samtools view -hb ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.sorted.bam -L test.bed > ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.bam
    samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.sorted.bam ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.bam
    samtools index ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.sorted.bam

    # python lrft2.py ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.sorted.bam ${OUT_PATH}/SMS ${ANNO_PATH}/${GENOME}/dm6_clean.transposon.2.mmi ${ANNO_PATH}/${GENOME}/dm6_clean.transposon.size 300


    python SMS.py ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.sorted.bam ${OUT_PATH}/SMS ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.mmi ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.size 300 sms.test ${ANNO_PATH}/${GENOME}/line_28_template.fa ${ANNO_PATH}/${GENOME}/line_28_template.mmi ${REPEATMASKER_FILE}
    # bedtools intersect -wa -wb -a ${OUT_PATH}/SMS/SMS.insertion.bed -b ${ANNO_PATH}/dm3/dm3.rmsk.bed | awk '{split($4,a,"|");for(i=1;i<=length(a);i++){split(a[i],b,":");if(b[1]==$13){print $0}}}' > ${OUT_PATH}/SMS/SMS.ex_rmsk.bed
    # bedtools intersect -v -a ${OUT_PATH}/SMS/SMS.insertion.bed -b ${OUT_PATH}/SMS/SMS.ex_rmsk.bed > ${OUT_PATH}/SMS/SMS.final.bed


    end=$(date +%s)
    time2=$(( $end - $start ))
    # echo -e "\n-----------\n" >> time_test_result.txt
    # echo -e $(date) >> time_test_result.txt
    # echo "${dep}" >> time_test_result.txt
    # echo ${time2}'s' >> time_test_result.txt
done



