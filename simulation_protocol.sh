#! /bin/bash


######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-d <dir>\tworking directory"
    echo -e "\t-n <name(s)>\tinbred line name(s) (id) OR sample name."
    echo -e "\t-r <fasta>\tindexed reference genome FASTA. (use this without -s)"
    echo -e "\t-s <F/M>\tsample sex, F/M."
    echo -e "\t-f <fasta>\tindexed female reference genome FASTA. (use this with -s)"
    echo -e "\t-m <fasta>\tindexed male reference genome FASTA. (use this with -s)"
    echo -e "\t-v <path>\tdirectory to indexed SNVs VCFs with genotype (GGVP for GRCh38; DGRP for dm3)."
    echo -e "\t-t <fasta>\tindexed TE CCS FASTA."
    echo -e "\t-N <int>\tnumber of insertions to generate for each TE."
    echo -e "\t-G <int>\tnumber length gradient."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":d:n:r:s:f:m:v:t:N:G:h" OPTION; do
    case $OPTION in
        d)  WORKING_DIR=$OPTARG;;
        n)  NAMEs=($OPTARG);;
        r)  FASTA=$OPTARG;;
        s)  SEX=$OPTARG;;
        f)  F_FASTA=$OPTARG;;
        m)  M_FASTA=$OPTARG;;
        v)  VCF_PATH=$OPTARG;;
        t)  TE_FASTA=$OPTARG;;
        N)  INS_NUM=$OPTARG;;
        G)  GRADIENT=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


### checking ###
echo -e "WORKING DIR:\t${WORKING_DIR}"
echo -e "NAME:\t${NAMEs[*]}"
echo -e "FASTA:\t${FASTA}"
echo -e "SEX:\t${SEX}"
echo -e "F_FASTA:\t${F_FASTA}"
echo -e "M_FASTA:\t${M_FASTA}"
echo -e "VCF_PATH:\t${VCF_PATH}"
echo -e "TE_FASTA:\t${TE_FASTA}"
echo -e "INS_NUM:\t${INS_NUM}"
echo -e "GRADIENT:\t${GRADIENT}"
### checking ###


### VARIABLES DECLARATION ###
source /data/tusers/zhongrenhu/Software/anaconda3/etc/profile.d/conda.sh


### SIMULATION WITH UN-PHASED SNV ###
if [ -z ${SEX} ];then
    echo -e "[ SIMULATION WITH UN-PHASED SNV ]"
    cd ${WORKING_DIR} && echo -e "[ NOW WORKING DIR:\t`pwd` ]"

    for NAME in ${NAMEs[*]}
    do
        ## BUILD TEMPLATE GENOME ##
        echo -e "[ BUILD TEMPLATE GENOME FOR ${NAME} ]"
        if [ -f ${NAME}_templateswithsnp/${NAME}_template.fa.fai ];then
            echo -e "[ TEMPLATE OF ${NAME} ALREADY EXSIST, SKIP. ]"
        else
            echo -e "[ CDM:\t/data/tusers/zhongrenhu/for_SMS/test/build_template_genome_unphased.sh -r ${FASTA} -v ${VCF_PATH} -n ${NAME} ]"
            conda activate bcftools_env
            /data/tusers/zhongrenhu/for_SMS/test/build_template_genome_unphased.sh -r ${FASTA} -v ${VCF_PATH} -n ${NAME}
        fi

        ## GENERATE FULL AND PARTIAL LENGTH INSERTION FROM TEMPLATE ##
        echo -e "[ GENERATE FULL AND PARTIAL LENGTH INSERTION FROM ${NAME} TEMPLATE ]"
        if [ -f ${NAME}_templateswithsnp/${NAME}.simulated_sv.summary ];then
            echo -e "[ INSERTIONS FOR ${NAME} ALREADY EXIST, SKIP. ]"
        else
            echo -e "[ GENERATE INSERTIONS FROM ${NAME}_templateswithsnp/${NAME}_template.fa START. ]"
            cd ${NAME}_templateswithsnp && echo -e "[ CMD:\tRscript /data/tusers/zhongrenhu/for_SMS/test/simulate_sv_genome.R ${TE_FASTA} ${NAME}_template.fa ${INS_NUM} ${GRADIENT} ${NAME} ]"
            conda activate R_env
            Rscript /data/tusers/zhongrenhu/for_SMS/test/simulate_sv_genome.R ${TE_FASTA} ${NAME}_template.fa ${INS_NUM} ${GRADIENT} ${NAME}
            cd ../
            echo -e "[ GENERATE INSERTIONS FROM ${NAME}_templateswithsnp/${NAME}_template.fa FINISH. ]"
        fi

        ## GENERATE SUBSETS OF INSERTIONS AND CALAULATE FREQUENCY ##
        echo -e "[ GENERATE SUBSETS OF INSERTIONS FROM ${NAME}_templateswithsnp/${NAME}.simulated_sv.summary ]"
        if [ -f ${NAME}_templateswithsnp/${NAME}.groundtruth.summary.seq ];then
            echo -e "[ SUBSETS OF ${NAME} ALREADY EXIST, SKIP. ]"
        else
            cd ${NAME}_templateswithsnp
            NROW=`wc -l ${NAME}.simulated_sv.summary`
            NUM=`awk -vN=${NROW} 'BEGIN{printf "%.0f", 0.5*N}'`
            /data/tusers/zhongrenhu/for_SMS/test/subset_ins.sh -n ${NAME} -M ${NUM} -R ${NAME}.simulated_sv.summary -S 4
            cd ..
        fi
    done
fi


### SIMULATION WITH PHASED SNV ###
# if [ ! -z ${SEX} ];then
#     echo -e "SIMULATION WITH PHASED SNV"
# fi


# if [ ! -f  ];then
# else
# fi