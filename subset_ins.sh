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


### checking ###
echo -e "NAME:\t${NAME}"
echo -e "N_INS:\t${N_INS}"
echo -e "REF_SUMMARY:\t${REF_SUMMARY}"
echo -e "SEX:\t${SEX}"
echo -e "G_SIZE:\t${G_SIZE}"
### checking ###


### GENERATE A GROUP OF SUBSET(S) AND CALCULATE INSERTION FREQUENCY FOR EACH GROUP ###
if [ ! -z ${SEX} ];then
    echo -e "[ GENERATE SUBSET(S) FOR ${NAME} START. ]"
    # GENERATE SUBSET(S)
    for ((j=1; j<=${G_SIZE}; j++))
    do
        echo -e "[ CMD:\tshuf -n ${N_INS} ${REF_SUMMARY} | sort -k1,1 -k2,2n > ./${NAME}.${j}.groundtruth.summary ]"
        shuf -n ${N_INS} ${REF_SUMMARY} | sort -k1,1 -k2,2n > ./${NAME}.${j}.groundtruth.summary
        echo -e "[ CMD:\tcut -f 1-6 ./${NAME}.${j}.groundtruth.summary > ./${NAME}.${j}.groundtruth.bed ]"
        cut -f 1-6 ./${NAME}.${j}.groundtruth.summary > ./${NAME}.${j}.groundtruth.bed
    done

    # CALCULATING FREQUENCY
    cat ./${NAME}.*.groundtruth.summary | sort -k1,1 -k2,2n | uniq -c > ./tmp.uniq
    echo -e "chrom\tstart\tend\tname\tadd_len\tstrand\tfrequency\tTE\tte_start\tte_end\ttsd_start\ttsd_end" > ./${NAME}.groundtruth.summary.bed
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,$9,$1,$10,$11,$12,$14,$15}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.bed
    echo -e "name\tte_seq\ttsd_seq\tins_seq" > ./${NAME}.groundtruth.summary.seq
    awk 'BEGIN{OFS="\t"} {print $8,$13,$16,$6}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.seq
    rm ./tmp.uniq
fi

# calculating frequency
cat ./${NAME}.*.groundtruth.summary | sort -k1,1 -k2,2n | uniq -c > ./tmp.uniq
echo -e "chrom\tstart\tend\tname\tadd_len\tstrand\tfrequency\tTE\tte_start\tte_len\ttsd_start\ttsd_len" > ./${NAME}.groundtruth.summary.bed
awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,$9,$1,$10,$11,$12,$14,$15}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.bed
echo -e "name\tte_seq\ttsd_seq\tins_seq" > ./${NAME}.groundtruth.summary.seq
awk 'BEGIN{OFS="\t"} {print $8,$13,$16,$6}' ./tmp.uniq >> ./${NAME}.groundtruth.summary.seq
rm ./tmp.uniq