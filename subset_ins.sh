#! /bin/bash

# usage: subset_ins.sh n_group g_size n_ins path_sv.summary out_path
# subset_ins.sh 2 5 1000 ./simulated_sv.summary ./

# geting parameters
n_group=$1
g_size=$2
n_ins=$3
ref_path=$4
out_path=$5

# parameter initialization
[ -z $n_group ] && n_group=1
[ -z $g_size ] && g_size=5
[ -z $n_ins ] && n_ins=1000
[ -z $ref_path ] && ref_path="/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/simulated_sv.summary"
[ -z $out_path ] && out_path="./"

# generating group(s) of subset(s) AND calculating insertion frequency for each group
for ((i=1; i<=$n_group; i++))
do
    out_dir=${out_path%*/}/group$i && mkdir $out_dir

    # generating subset(s)
    for ((j=1; j<=$2; j++))
    do
        shuf -n $n_ins $ref_path | sort -k1,1 -k2,2n > $out_dir/group$i-sub$j.summary
        cut -f 1-6 $out_dir/group$i-sub$j.summary > $out_dir/group$i-sub$j.bed
    done

    # calculating frequency
    cat $out_dir/group$i-sub*summary | sort -k1,1 -k2,2n | uniq -c > $out_dir/tmp.uniq
    echo -e "chrom\tstart\tend\tname\tadd_len\tstrand\tfrequency\tTE\tte_start\tte_len\ttsd_start\ttsd_len" > $out_dir/group$i-ins-summary.bed
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,$9,$1,$10,$11,$12,$14,$15}' $out_dir/tmp.uniq >> $out_dir/group$i-ins-summary.bed
    echo -e "name\tte_seq\ttsd_seq\tins_seq" > $out_dir/group$i-ins-seq
    awk 'BEGIN{OFS="\t"} {print $8,$13,$16,$6}' $out_dir/tmp.uniq >> $out_dir/group$i-ins-seq
    rm $out_dir/tmp.uniq

done