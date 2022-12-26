#!/bin/bash
start=$(date +%s)

# DATA_PATH="/data/tusers/boxu/lrft2/rawdata/simu_div0/downsample"
DATA_PATH="/data/tusers/boxu/lrft2/rawdata/simu_div0/downsample_by_reads/"
# OUT_PATH_dir="/data/tusers/boxu/lrft2/result/"
OUT_PATH_dir="/data/tusers/boxu/lrft2/validation/"
lrft_BIN_PATH="/data/tusers/boxu/lrft/scripts"
ANNO_PATH="/data/tusers/boxu/annotation"
GENOME="dm6_clean"
# SAMPLE="line_28.genome"
SAMPLE="line_28"
CPU="10"
meth="lrft"

BWA_INDEX="/data/tusers/boxu/annotation/dm3/BWAIndex"


TRANSPOSON_INDEX="/data/tusers/boxu/annotation/dm3/dm3.transposon_for_simulaTE.mmi"
TRANSPOSON_SIZE="/data/tusers/boxu/annotation/dm3/dm3.transposon_for_simulaTE.size"
GENOME_FA="/data/tusers/boxu/annotation/dm3/line_28_template.fa"
GENOME_INDEX="/data/tusers/boxu/annotation/dm3/line_28_template.mmi"
TE_ANNO_FA="/data/tusers/boxu/annotation/dm3/dm3.transposon_for_simulaTE.fa"
REPEATMASKER_FILE="/data/tusers/boxu/annotation/dm3/dm3.rmsk.bed"
# python /data/tusers.ds/boxu/lrft2/scripts/setup.py build && cp /data/tusers.ds/boxu/lrft2/scripts/build/lib.linux-x86_64-3.7/read_alignment.cpython-37m-x86_64-linux-gnu.so /data/tusers.ds/boxu/lrft2/bin/
# chr2L:7,501,451-7,502,704
# OUT_PATH=${OUT_PATH_dir}
[ ! -d ${OUT_PATH_dir} ] && mkdir ${OUT_PATH_dir}

for dep in 50 # 1 2 3 4 5 10 20 30 40 50
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
    
    [ ! -d ${OUT_PATH}/SMS2 ] && mkdir ${OUT_PATH}/SMS2
    [ ! -d ${OUT_PATH}/SMS2/temp ] && mkdir ${OUT_PATH}/SMS2/temp

    [ ! -d ${OUT_PATH}/TEMP2 ] && mkdir ${OUT_PATH}/TEMP2
    [ ! -d ${OUT_PATH}/tldr ] && mkdir ${OUT_PATH}/tldr

    
    # NGS mapping
    if [ ! -f ${DATA_PATH}/${SAMPLE}_ngs_${dep}X.sort.bam ];then
        LEFT=${DATA_PATH}/${SAMPLE}_${dep}X_1.fastq
        RIGHT=${DATA_PATH}/${SAMPLE}_${dep}X_2.fastq

        bwa mem -t ${CPU} ${BWA_INDEX}/line_28_template ${LEFT} ${RIGHT} > ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sam 2>${DATA_PATH}/${SAMPLE}_ngs.${dep}X.bwamem.log
        samtools view -bhS -@ ${CPU} ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sam > ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.bam
        samtools sort -@ ${CPU} -o ${DATA_PATH}/${SAMPLE}_ngs_${dep}X.sort.bam ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.bam
        samtools index -@ ${CPU} ${DATA_PATH}/${SAMPLE}_ngs_${dep}X.sort.bam
    fi

    if [ ! -f ${OUT_PATH}/${SAMPLE}_ngs_${dep}X_q0.sorted.bam ];then
        samtools view -bhS -q 5 ${DATA_PATH}/${SAMPLE}_ngs_${dep}X.sort.bam > ${OUT_PATH}/${SAMPLE}_ngs_sort_${dep}X_q0.bam
        samtools sort -@ ${CPU} -o ${OUT_PATH}/${SAMPLE}_ngs_${dep}X_q0.sorted.bam ${OUT_PATH}/${SAMPLE}_ngs_sort_${dep}X_q0.bam
        samtools index -@ ${CPU} ${OUT_PATH}/${SAMPLE}_ngs_${dep}X_q0.sorted.bam
    fi



    # TGS mapping
    if [ ! -f ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.bam ];then
        minimap2 -aYx map-pb ${ANNO_PATH}/dm3/line_28_template.mmi ${DATA_PATH}/${SAMPLE}_pacbio_${dep}X.fasta > ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.sam
        samtools sort -@ 8 -O bam -o  ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.bam  ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.sam
        samtools index -@ 8 ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.bam
    fi


    if [ ! -f ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam ];then
        samtools view -bhSq 5 ${DATA_PATH}/${SAMPLE}_tgs.${dep}X.bam > ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.bam
        samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.bam
        samtools index -@ 10 ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam
    fi
    
    
    # if [ ! -f ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.sorted.bam ];then
    #     samtools view -bhSq 5 ${DATA_PATH}/${SAMPLE}.genome.${dep}X.bam > ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.bam
    #     samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.sorted.bam ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.bam
    #     samtools index -@ 10 ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.sorted.bam
    # fi



    # NGS 
    # TEMP2
    #echo "TEMP2 insertion -l ${DATA_PATH}/${SAMPLE}_${dep}X_1.fastq -r ${DATA_PATH}/${SAMPLE}_${dep}X_2.fastq -i ${OUT_PATH}/${SAMPLE}_ngs.${dep}X.q0.sorted.bam -g ${ANNO_PATH}/dm3/line_28_template.fa -R ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.fa -t ${ANNO_PATH}/dm3/dm3.rmsk.bed -o ${OUT_PATH}/TEMP2 -p ${SAMPLE}_ngs.${dep}X -d -c ${CPU}"
    source /home/boxu/bin/anaconda3/bin/activate py2 && TEMP2 insertion -l ${DATA_PATH}/${SAMPLE}_${dep}X_1.fastq -r ${DATA_PATH}/${SAMPLE}_${dep}X_2.fastq -i ${DATA_PATH}/${SAMPLE}_ngs_${dep}X.sort.bam -g ${ANNO_PATH}/dm3/line_28_template.fa -R ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.fa -t ${ANNO_PATH}/dm3/dm3.rmsk.bed -o ${OUT_PATH}/TEMP2 -p ${SAMPLE}_ngs.${dep}X -d -c ${CPU} && conda deactivate

    # tldr
    source /home/boxu/bin/anaconda3/bin/activate tldr && tldr -m 1 --max_cluster_size 600 --max_te_len 12000 --min_te_len 150 --embed_minreads 0 -b ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam -e ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.2.fa -r ${ANNO_PATH}/dm3/line_28_template.fa --color_consensus -o ${OUT_PATH}/tldr/${SAMPLE}_pacbio && conda deactivate
    

    # SMS
    # echo -e "chr2L\t13851548\t13868000" > test.bed
    # samtools view -hb ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.sorted.bam -L test.bed > ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.region.bam
    # samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.region.sorted.bam ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.region.bam
    # samtools index ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.region.sorted.bam

    # python lrft2.py ${OUT_PATH}/${SAMPLE}.genome.${dep}X.q0.region.sorted.bam ${OUT_PATH}/SMS ${TRANSPOSON_INDEX} ${TRANSPOSON_SIZE} 1500 SMS.repeat.test4 ${GENOME_FA} ${GENOME_INDEX} ${REPEATMASKER_FILE}

    prefix_re="SMS"
    # lrft2.py = SMS.py

    # python lrft2.py ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam ${OUT_PATH}/SMS ${TRANSPOSON_INDEX} ${TRANSPOSON_SIZE} 1500 ${prefix_re} ${GENOME_FA} ${GENOME_INDEX} ${TE_ANNO_FA} ${REPEATMASKER_FILE}
    python lrft2.py -b ${OUT_PATH}/${SAMPLE}_tgs.${dep}X.q0.sorted.bam -o ${OUT_PATH}/SMS -t ${TRANSPOSON_INDEX} -s ${TRANSPOSON_SIZE} -f 1500 -p ${prefix_re} -g ${GENOME_FA} -i ${GENOME_INDEX} -a ${TE_ANNO_FA} -r ${REPEATMASKER_FILE}

    # grep pass ${OUT_PATH}/SMS/${prefix_re}.insertion.bed | grep -v nest | grep -v TART | grep -v del | awk '{split($14,a,"_");split($8,b,"_");dep=a[1]+a[2]+a[3];print b[2],b[3],$9,dep,$15}' > ${OUT_PATH}/SMS/${prefix_re}.${dep}.txt
    # mv ${OUT_PATH}/SMS/${prefix_re}.${dep}.txt /home/boxu/temp/wtdbg2/for_qc/2/

    # python lrft2.py ${OUT_PATH}/${SAMPLE}.${dep}X.q0.sorted.bam ${OUT_PATH}/SMS ${ANNO_PATH}/${GENOME}/dm6_clean.transposon.2.mmi ${ANNO_PATH}/${GENOME}/dm6_clean.transposon.size 300
    # bedtools intersect -wa -wb -a ${OUT_PATH}/SMS/SMS.insertion.bed -b ${ANNO_PATH}/dm3/dm3.rmsk.bed | awk '{split($4,a,"|");for(i=1;i<=length(a);i++){split(a[i],b,":");if(b[1]==$13){print $0}}}' > ${OUT_PATH}/SMS/SMS.ex_rmsk.bed
    # bedtools intersect -v -a ${OUT_PATH}/SMS/SMS.insertion.bed -b ${OUT_PATH}/SMS/SMS.ex_rmsk.bed > ${OUT_PATH}/SMS/SMS.final.bed


    end=$(date +%s)
    time2=$(( $end - $start ))
    # echo -e "\n-----------\n" >> time_test_result.20X.txt
    # echo -e $(date) >> time_test_result.20X.txt
    # echo "${dep}" >> time_test_result.20X.txt
    # echo ${time2}'s' >> time_test_result.20X.txt
done



exit

if [ ! -f ${OUT_PATH}/${SAMPLE}.genome.bam ];then

    minimap2  -t 10 -aYx map-ont ${ANNO_PATH}/${GENOME}/line_28_template.mmi ${DATA_PATH}/${SAMPLE}.fasta > ${OUT_PATH}/${SAMPLE}.genome.sam
    samtools sort -@ 8 -O bam -o ${OUT_PATH}/${SAMPLE}.genome.bam ${OUT_PATH}/${SAMPLE}.genome.sam
    samtools index ${OUT_PATH}/${SAMPLE}.genome.bam
fi

if [ ! -f ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ];then
    samtools view -bhSq 0 ${OUT_PATH}/${SAMPLE}.genome.bam > ${OUT_PATH}/${SAMPLE}.q0.bam
    samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ${OUT_PATH}/${SAMPLE}.q0.bam
    samtools index -@ 10 ${OUT_PATH}/${SAMPLE}.q0.sorted.bam
fi


# samtools view -bhSq 0 ${OUT_PATH}/${SAMPLE}.q0.sorted.bam -L test.bed > ${OUT_PATH}/${SAMPLE}.q0.region.bam
# samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}.q0.region.sorted.bam ${OUT_PATH}/${SAMPLE}.q0.region.bam
# samtools index ${OUT_PATH}/${SAMPLE}.q0.region.sorted.bam

# python lrft2.py ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ${OUT_PATH}/SMS ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.mmi ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.size 300

# python lrft2.py ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ${OUT_PATH}/SMS ${TRANSPOSON_INDEX} ${TRANSPOSON_SIZE} 300 SMS ${GENOME_FA} ${GENOME_INDEX} ${REPEATMASKER_FILE}

echo -e  "chr2L\t15938978\t15938989" > test.bed
samtools view -hb ${OUT_PATH}/${SAMPLE}.q0.sorted.bam -L test.bed > ${OUT_PATH}/${SAMPLE}.q0.region.bam
samtools sort -@ 10 -o ${OUT_PATH}/${SAMPLE}.q0.region.sorted.bam ${OUT_PATH}/${SAMPLE}.q0.region.bam
samtools index ${OUT_PATH}/${SAMPLE}.q0.region.sorted.bam

python lrft2.py ${OUT_PATH}/${SAMPLE}.q0.region.sorted.bam ${OUT_PATH}/SMS ${TRANSPOSON_INDEX} ${TRANSPOSON_SIZE} 300 SMS.test ${GENOME_FA} ${GENOME_INDEX} ${REPEATMASKER_FILE}


    #echo -e "chr2R\t1037043\t1075394" > test.bed
    #echo -e "chrU\t3204103\t3208590" > test.bed # 这段比较重要吧，空窗期
    #echo -e "chr2L\t22434048\t22437297" > test.bed
    # echo -e "chr2L\t1684519\t1692744" > test.bed # l_black_num == black_num
    # echo -e "chr2L\t7792778\t7795650" > test.bed
    #echo -e "chr2L\t19547988\t19548305" > test.bed
    #echo -e "chrU\t2971477\t3086434" > test.bed
    #echo -e "chr2RHet\t60175\t60183" > test.bed # 没有double clip read，只有一端clip 的supporting reads
    #echo -e "chrU\t4001230\t4058708" >> test.bed
    # echo -e "chrU\t8140745\t8140947" > test.bed
    # echo -e "chr2R\t14773963\t14782414" > test.bed
    # # echo -e "chr3L\t20814897\t20819283" > test.bed
    # echo -e "chr3L\t23205501\t23205847" > test.bed
