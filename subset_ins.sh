#! /bin/bash
# usage: subset_ins.sh n_group g_size n_ins path_sv.summary out_path
# subset_ins.sh 2 5 1000 ./simulated_sv.summary ./

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-n <name>\tinbred line name(id) OR sample name."
    echo -e "\t-M <int>\tinertions number of each subset that use to generate target haplotype."
    echo -e "\t-R <ref.summary>\ttotal simulated insertions that used for generate subsets."
    echo -e "\t-s <F/M>\tsample sex, F/M."
    echo -e "\t-S <int>\tgroup size."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":n:M:R:s:S:h" OPTION; do
    case $OPTION in
        n)  NAME=$OPTARG;;
        M)  N_INS=$OPTARG;;
        R)  REF_SUMMARY=$OPTARG;;
        s)  SEX=$OPTARG;;
        S)  G_SIZE=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

export PATH=$PATH:/data/tusers/zhongrenhu/Software/anaconda3/bin/


### Checking ###
echo -e "NAME:\t${NAME}"
echo -e "N_INS:\t${N_INS}"
echo -e "REF_SUMMARY:\t${REF_SUMMARY}"
echo -e "SEX:\t${SEX}"
echo -e "G_SIZE:\t${G_SIZE}"
### CHEKING ###


### GENERATE A GROUP OF SUBSET(S) AND CALCULATE INSERTION FREQUENCY FOR EACH GROUP ###
if [ -z ${SEX} ];then
    echo -e "[ GENERATE SUBSET(S) FOR ${NAME} START. ]"

    for ((j=1; j<=${G_SIZE}; j++))
    do
        # GENERATE SUBSET(S)
        if [ ! -f ${NAME}.homozygous.groundtruth.summary ];then
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H
            awk '$5 !~ /[ATCG]*[SN]+[ATCG]*/{print $0}' ${NAME}.homozygous.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.homozygous.groundtruth.summary
            awk '$5 !~ /[ATCG]*[SN]+[ATCG]*/{print $0}' ${NAME}.${j}.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
            cut -f 1-6 ${NAME}.homozygous.groundtruth.summary > ${NAME}.homozygous.groundtruth.bed
        else
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY}
            awk '$5 !~ /[ATCG]*[SN]+[ATCG]*/{print $0}' ${NAME}.${j}.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
        fi

        cut -f 1-6 ${NAME}.${j}.groundtruth.summary > ${NAME}.${j}.groundtruth.bed
        cat ${NAME}.homozygous.groundtruth.summary >> tmp.all
        cat ${NAME}.${j}.groundtruth.summary >> tmp.all

        # INTEGRATE INSERTIONS INTO TEMPLATE GENOME
        if [ ! -f ${NAME}.${j}.fa.fai ];then
            echo -e "[ INTEGRATE ${NAME}.homozygous.groundtruth.bed AND ${NAME}.${j}.groundtruth.bed INTO ${NAME}_template.fa START. ]"
            mkdir tmp && cat ${NAME}.homozygous.groundtruth.bed ${NAME}.${j}.groundtruth.bed > tmp.bed
            VISOR HACk -g ${NAME}_template.fa -b tmp.bed -o tmp
            mv tmp/h1.fa ./${NAME}.${j}.fa && mv tmp/h1.fa.fai ./${NAME}.${j}.fa.fai
            rm -r tmp && rm tmp.bed
            echo -e "[ INTEGRATION FOR ${NAME}.${j} FINISH. ]"
        fi
    done

    # CALCULATING FREQUENCY
    sort -k1,1 -k2,2n ./tmp.all | uniq -c > ./tmp.uniq
    echo -e "chrom\tstart\tend\tname\tadd_len\tstrand\tfrequency\tTE\tte_start\tte_end\ttsd_start\ttsd_end\tn_mutates\tFullLength" > ./${NAME}.groundtruth.summary.bed
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,$9,$1,$10,$11,$12,$14,$15,$17,$19}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.bed
    echo -e "name\tte_seq\ttsd_seq\tins_seq\torigin_seq" > ./${NAME}.groundtruth.summary.seq
    awk 'BEGIN{OFS="\t"} {print $8,$13,$16,$6,$18}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.seq
    rm tmp.uniq && rm tmp.all
    echo -e "[ GENERATE SUBSET(S) FOR ${NAME} FINISH. ]"
fi