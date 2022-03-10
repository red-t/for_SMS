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
        if [ ! -f ${NAME}.homozygous.groundtruth.bed ];then
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H
            awk '$5 ~ /^[ATCG]*$/{print $0}' ${NAME}.homozygous.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.homozygous.groundtruth.summary
            awk '$5 ~ /^[ATCG]*$/{print $0}' ${NAME}.${j}.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
            cut -f 1-6 ${NAME}.homozygous.groundtruth.summary > ${NAME}.homozygous.groundtruth.bed
            cut -f 1-6 ${NAME}.${j}.groundtruth.summary > ${NAME}.${j}.groundtruth.bed
        elif [ ! -f ${NAME}.${j}.groundtruth.bed ];then
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY}
            awk '$5 ~ /^[ATCG]*$/{print $0}' ${NAME}.${j}.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
            cut -f 1-6 ${NAME}.${j}.groundtruth.summary > ${NAME}.${j}.groundtruth.bed
        fi

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


if [ ! -z ${SEX} ];then
    echo -e "[ GENERATE SUBSET(S) FOR ${NAME} START. ]"

    for ((j=1; j<=${G_SIZE}; j++))
    do
        # GENERATE SUBSET(S)
        if [ ! -f ${NAME}.homozygous.groundtruth.h2.bed ];then
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} -H
            awk '$5~/^[ATCG]*$/ && $19~/^[ATCG]*$/{print $0}' ${NAME}.homozygous.groundtruth.summary | sort -k1,1 -k2,2n > tmp.gt && mv tmp.gt ${NAME}.homozygous.groundtruth.summary
            awk '$5~/^[ATCG]*$/ && $19~/^[ATCG]*$/{print $0}' ${NAME}.${j}.groundtruth.summary > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' ${NAME}.homozygous.groundtruth.summary > ${NAME}.homozygous.groundtruth.h1.bed
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$12$19,$6}' ${NAME}.homozygous.groundtruth.summary > ${NAME}.homozygous.groundtruth.h2.bed

            NR=`cat ${NAME}.${j}.groundtruth.summary | wc -l`; NR1=`awk -vN=${NR} 'BEGIN{printf "%.0f",0.5*N}'`; NR2=$((NR-NR1))
            head -n ${NR1} ${NAME}.${j}.groundtruth.summary | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > ${NAME}.${j}.groundtruth.h1.bed
            tail -n ${NR2} ${NAME}.${j}.groundtruth.summary | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$12$19,$6}' > ${NAME}.${j}.groundtruth.h2.bed
        elif [ ! -f ${NAME}.${j}.groundtruth.bed ];then
            echo -e "[ CMD:\tpython /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY} ]"
            python /data/tusers/zhongrenhu/for_SMS/test/my_shuf.py -p ${NAME}.${j} -M ${N_INS} -R ${REF_SUMMARY}
            awk '$5~/^[ATCG]*$/ && $19~/^[ATCG]*$/{print $0}' ${NAME}.${j}.groundtruth.summary > tmp.gt && mv tmp.gt ${NAME}.${j}.groundtruth.summary
            
            NR=`cat ${NAME}.${j}.groundtruth.summary | wc -l`; NR1=`awk -vN=${NR} 'BEGIN{printf "%.0f",0.5*N}'`; NR2=$((NR-NR1))
            head -n ${NR1} ${NAME}.${j}.groundtruth.summary | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > ${NAME}.${j}.groundtruth.h1.bed
            tail -n ${NR2} ${NAME}.${j}.groundtruth.summary | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$12$19,$6}' > ${NAME}.${j}.groundtruth.h2.bed
        fi

        cat ${NAME}.homozygous.groundtruth.summary >> tmp.all && cat ${NAME}.homozygous.groundtruth.summary >> tmp.all
        cat ${NAME}.${j}.groundtruth.summary >> tmp.all

        # INTEGRATE INSERTIONS INTO TEMPLATE GENOME
        if [ ! -f ${NAME}.${j}.h2.fa.fai ];then
            echo -e "[ INTEGRATION FOR ${NAME}.${j} START. ]"
            # INTEGRATION FOR H1
            mkdir tmp
            cat ${NAME}.homozygous.groundtruth.h1.bed ${NAME}.${j}.groundtruth.h1.bed > tmp.bed
            VISOR HACk -g ${NAME}_template_h1.fa -b tmp.bed -o tmp
            mv tmp/h1.fa ./${NAME}.${j}.h1.fa && mv tmp/h1.fa.fai ./${NAME}.${j}.h1.fa.fai

            # INTEGRATION FOR H2
            rm tmp.bed
            cat ${NAME}.homozygous.groundtruth.h2.bed ${NAME}.${j}.groundtruth.h2.bed > tmp.bed
            VISOR HACk -g ${NAME}_template_h2.fa -b tmp.bed -o tmp
            mv tmp/h1.fa ./${NAME}.${j}.h2.fa && mv tmp/h1.fa.fai ./${NAME}.${j}.h2.fa.fai
            rm -r tmp && rm tmp.bed
            echo -e "[ INTEGRATION FOR ${NAME}.${j} FINISH. ]"
        fi
    done

    # CALCULATING FREQUENCY
    sort -k1,1 -k2,2n ./tmp.all | uniq -c > ./tmp.uniq
    echo -e "chrom\tstart\tend\tname\tadd_len\tstrand\tfrequency\tTE\tte_start\tte_end\ttsd_start\ttsd_end\tn_mutates\tFullLength" > ./${NAME}.groundtruth.summary.bed
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,$9,$1,$10,$11,$12,$14,$15,$17,$19}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.bed
    echo -e "name\tte_seq\ttsd_seq\th2_tsd_seq\tins_seq\torigin_seq" > ./${NAME}.groundtruth.summary.seq
    awk 'BEGIN{OFS="\t"} {print $8,$13,$16,$20,$6,$18}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.seq
    rm tmp.uniq && rm tmp.all
    echo -e "[ GENERATE SUBSET(S) FOR ${NAME} FINISH. ]"
fi