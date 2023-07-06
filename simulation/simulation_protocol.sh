#! /bin/bash

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-d <dir>\tworking directory"
    echo -e "\t-r <fasta>\tindexed reference genome FASTA."
    echo -e "\t-t <fasta>\tindexed TE CCS FASTA."
    echo -e "\t-N <int>\tnumber of genome to simulate."
    echo -e "\t-R <float>\tTE sequence divergence rate (percent)."
    echo -e "\t--sub-N <int>\tsub pop size, number of genome that used to generate reads."
    echo -e "\t--germline-count <int>\tinsertions numbers of each frequency gradients."
    echo -e "\t--avg-somatic-count <int>\taverage somatic insertion mumber in each genome."
    echo -e "\t--min-distance <int>\tminimum distance between TE insertions."
    echo -e "\t--depth <int>\tsequencing depth of the NGS/TGS data."
    echo -e "\t--ngs-len <int>\tread length of NGS reads."
    echo -e "\t--ngs-inner <int>\tinner size of NGS fragment."
    echo -e "\t--ngs-std <int>\tStanard derivation of NGS inner size length distribution."
    echo -e "\t--ngs-err <float>\terror rate of NGS reads (fraction)."
    echo -e "\t--tgs-maxl <int>\tMax length of TGS reads."
    echo -e "\t--tgs-minl <int>\tMin length of TGS reads."
    echo -e "\t--tgs-len <int>\tMean length of TGS reads."
    echo -e "\t--tgs-alpha <float>\talpha of gamma distribution model applied for TGS length generation."
    echo -e "\t--tgs-loc <float>\tloc of gamma distribution model applied for TGS length generation."
    echo -e "\t--tgs-beta <float>\tbeta of gamma distribution model applied for TGS length generation."
    echo -e "\t--tgs-err <float>\terror rate of TGS reads (fraction)."
    echo -e "\t--err-frac \tthe fraction of each type of errors, mis:del:ins."
    echo -e "\t-h \tShow this information"
}


######## Getting Parameters ########
ARGS=`getopt -o d:r:t:N:R:h --long sub-N:,germline-count:,avg-somatic-count:,min-distance:,depth:,ngs-len:,ngs-inner:,ngs-std:,ngs-err:,tgs-maxl:,tgs-minl:,tgs-len:,tgs-alpha:,tgs-loc:,tgs-beta:,tgs-err:,err-frac: -n "$0" -- "$@"`
if [ $? != 0 ]; then
    echo "Terminating..."
    exit 1
fi

eval set -- "${ARGS}"

while true
do
    case "$1" in
        -d) 
            WORKING_DIR=$2
            shift 2
            ;;
        -r)
            REF_FA=$2
            shift 2
            ;;
        -t)
            TE_FA=$2
            shift 2
            ;;
        -N)
            POP_SIZE=$2
            shift 2
            ;;
        -R)
            D_RATE=$2
            shift 2
            ;;
        --sub-N)
            SUB_POP_SIZE=$2
            shift 2
            ;;
        --germline-count)
            G_COUNT=$2
            shift 2
            ;;
        --avg-somatic-count)
            AVG_S_COUNT=$2
            shift 2
            ;;
        --min-distance)
            MIN_DIST=$2
            shift 2
            ;;
        --depth)
            DEPTH=$2
            shift 2
            ;;
        --ngs-len)
            NGS_LEN=$2
            shift 2
            ;;
        --ngs-inner)
            NGS_INNER=$2
            shift 2
            ;;
        --ngs-std)
            NGS_STD=$2
            shift 2
            ;;
        --ngs-err)
            NGS_ERR=$2
            shift 2
            ;;
        --tgs-maxl)
            TGS_MAXL=$2
            shift 2
            ;;
        --tgs-minl)
            TGS_MINL=$2
            shift 2
            ;;
        --tgs-len)
            TGS_MEANL=$2
            shift 2
            ;;
        --tgs-alpha)
            TGS_ALPHA=$2
            shift 2
            ;;
        --tgs-loc)
            TGS_LOC=$2
            shift 2
            ;;
        --tgs-beta)
            TGS_BETA=$2
            shift 2
            ;;
        --tgs-err)
            TGS_ERR=$2
            shift 2
            ;;
        --err-frac)
            ERR_FRAC=$2
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -h)
            help_info
            shift
            exit 1
            ;;
        *)
            echo "Internal error!"
            help_info
            exit 1
            ;;
    esac
done


### Default parameters ###
[ -z $D_RATE ] && D_RATE=0
[ -z $NGS_LEN ] && NGS_LEN=150
[ -z $NGS_INNER ] && NGS_INNER=200
[ -z $NGS_STD ] && NGS_STD=20
[ -z $NGS_ERR] && NGS_ERR=0.0005
[ -z $TGS_MAXL] && TGS_MAXL=16000
[ -z $TGS_MINL] && TGS_MINL=6000
[ -z $TGS_MEANL] && TGS_MEANL=8690
[ -z $TGS_ALPHA] && TGS_ALPHA=3.3567764489878322
[ -z $TGS_LOC] && TGS_LOC=5999.934752151781
[ -z $TGS_BETA] && TGS_BETA=801.77408988051
[ -z $TGS_ERR] && TGS_ERR=0.1
[ -z $ERR_FRAC ] && ERR_FRAC='4:15:81'

### Checking Parameters ###
echo -e "simulation_protocol.sh parameters:"
echo -e "WORKING DIR:\t${WORKING_DIR}"
echo -e "REF_FA:\t${REF_FA}"
echo -e "TE_FA:\t${TE_FA}"
echo -e "POP_SIZE:\t${POP_SIZE}"
echo -e "D_RATE:\t${D_RATE}"
echo -e "SUB_POP_SIZE:\t${SUB_POP_SIZE}"
echo -e "G_COUNT:\t${G_COUNT}"
echo -e "AVG_S_COUNT:\t${AVG_S_COUNT}"
echo -e "MIN_DIST:\t${MIN_DIST}"
echo -e "DEPTH:\t${DEPTH}"
echo -e "NGS_LEN:\t${NGS_LEN}"
echo -e "NGS_INNER:\t${NGS_INNER}"
echo -e "NGS_STD:\t${NGS_STD}"
echo -e "NGS_ERR:\t${NGS_ERR}"
echo -e "TGS_MAXL:\t${TGS_MAXL}"
echo -e "TGS_MINL:\t${TGS_MINL}"
echo -e "TGS_MEANL:\t${TGS_MEANL}"
echo -e "TGS_ALPHA:\t${TGS_ALPHA}"
echo -e "TGS_LOC:\t${TGS_LOC}"
echo -e "TGS_BETA:\t${TGS_BETA}"
echo -e "TGS_ERR:\t${TGS_ERR}"
echo -e "ERR_FRAC:\t${ERR_FRAC}"
### Checking Parameters ###


### Variables Declaration ###
# source /data/tusers/zhongrenhu/Software/anaconda3/etc/profile.d/conda.sh
prog_path=`readlink -f $0` && prog_path=`dirname $prog_path`



### Generate Population Genome Definition (pgd-file) ###
S_COUNT=$(($AVG_S_COUNT * $POP_SIZE))
genome_size=`cut -f 2 ${REF_FA}.fai | awk 'BEGIN{n=0} {n=n+$1} END{print n}'`
contigs=(`cut -f 1 ${REF_FA}.fai`)
contigs_size=(`cut -f 2 ${REF_FA}.fai`)
contigs_count=`wc -l ${REF_FA}.fai | awk '{print $1}'`
del_g=0
del_s=0
for((i=0; i<$contigs_count; i++))
do
    # Calculate insertions numbers of each contig based on their length
    if [ $i -lt $(($contigs_count - 1)) ]; then
        tmp_g_count=`awk -v g_count=${G_COUNT} -v g_l=${genome_size} -v c_l=${contigs_size[$i]} 'BEGIN{ratio=c_l/g_l; print int(g_count*ratio+0.5)}'`
        tmp_s_count=`awk -v s_count=${S_COUNT} -v g_l=${genome_size} -v c_l=${contigs_size[$i]} 'BEGIN{ratio=c_l/g_l; print int(s_count*ratio+0.5)}'`
        del_g=$(($del_g + $tmp_g_count))
        del_s=$(($del_s + $tmp_s_count))
    elif [ $i -eq $(($contigs_count - 1)) ]; then
        tmp_g_count=$(($G_COUNT - $del_g))
        tmp_s_count=$(($S_COUNT - $del_s))
    fi
    
    # Extract contig sequence AND generate header & body of the pgd-file
    echo -e "Defining TE landscapes for ${contigs[$i]}:"
    if [ ! -f ${contigs[$i]}/${contigs[$i]}.ins.summary ]; then
        samtools faidx ${REF_FA} ${contigs[$i]} > ${contigs[$i]}.tmp.chasis.fasta && samtools faidx ${contigs[$i]}.tmp.chasis.fasta
        python $prog_path/my-define-landscape_random-insertions-freq-range.py --chassis ${contigs[$i]}.tmp.chasis.fasta --te-seqs ${TE_FA} --N ${POP_SIZE} --divergence-rate ${D_RATE} --germline-count ${tmp_g_count} --somatic-count ${tmp_s_count} --min-distance ${MIN_DIST}
        
        # Divide the pgd-file into several subsets
        if [ -f ${contigs[$i]}.tmp.pgd.header ]; then
            N_SUB=`awk -v sub_pop=$SUB_POP_SIZE -v pop=$POP_SIZE 'BEGIN{n_sub=int(pop/sub_pop); print n_sub}'`
            for((j=0; j<$N_SUB; j++))
            do
                start=`awk -v idx=$j -v pop=$SUB_POP_SIZE 'BEGIN{print idx*pop+2}'`
                end=`awk -v idx=$j -v pop=$SUB_POP_SIZE 'BEGIN{print (idx+1)*pop+1}'`
                cat ${contigs[$i]}.tmp.pgd.header > ${contigs[$i]}.$j.pgd
                cut -d " " -f 1,$start-$end ${contigs[$i]}.tmp.pgd.body >> ${contigs[$i]}.$j.pgd
            done
            rm ${contigs[$i]}.tmp.pgd.header && rm ${contigs[$i]}.tmp.pgd.body
        fi

        if [ -f ${contigs[$i]}.0.pgd ]; then
            mkdir ${contigs[$i]} && mv ${contigs[$i]}.* ${contigs[$i]}
        fi
    else
        echo -e "The population-genome-definition (pgd) files of ${contigs[$i]} already exist, skip."
    fi
done


### Build Population Genome & Generate 50X NGS, 3GS data ###
for((i=0; i<$contigs_count; i++))
do
    if [ -f ${contigs[$i]}/${contigs[$i]}.0.pgd ]; then
        cd ${contigs[$i]}
        N_SUB=`awk -v sub_pop=$SUB_POP_SIZE -v pop=$POP_SIZE 'BEGIN{n_sub=int(pop/sub_pop); print n_sub}'`
        NGS_READS=`awk -v depth=$DEPTH -v g_l=$genome_size -v c_l=${contigs_size[$i]} -v n_sub=$N_SUB -v r_l=150 'BEGIN{ratio=c_l/g_l; ngs_reads=int(ratio*(depth*g_l/(2*r_l))/n_sub); print ngs_reads}'`
        TGS_READS=`awk -v depth=$DEPTH -v g_l=$genome_size -v c_l=${contigs_size[$i]} -v n_sub=$N_SUB -v r_l=$TGS_MEANL 'BEGIN{ratio=c_l/g_l; tgs_reads=int(ratio*(depth*g_l/r_l)/n_sub); print tgs_reads}'`
        for((j=0; j<$N_SUB; j++))
        do
            # build population geneome
            echo -e "[ Build ${j}th sub population genome of ${contigs[$i]} ]"
            if [ ! -f ${contigs[$i]}.$j.fa ]; then
                # python $prog_path/simulate/build-population-genome.py --pgd ${contigs[$i]}.$j.pgd --te-seqs $TE_FA --chassis ${contigs[$i]}.tmp.chasis.fasta --output ${contigs[$i]}.$j.fa --ins-seq ${contigs[$i]}.ins.sequence --sub_idx $j --sub_size $SUB_POP_SIZE
                python $prog_path/build-population-genome.py --pgd ${contigs[$i]}.$j.pgd --te-seqs $TE_FA --chassis ${contigs[$i]}.tmp.chasis.fasta --output ${contigs[$i]}.$j.fa --ins-seq "" --sub_idx $j --sub_size $SUB_POP_SIZE
            fi
            # generate 50X NGS
            echo -e "[ Generate NGS data from ${j}th sub population genome of ${contigs[$i]} ]"
            if [ ! -f ${contigs[$i]}.${j}_2.fastq ]; then
                python $prog_path/read_pool-seq_illumina-PE.py --pg ${contigs[$i]}.$j.fa --read-length $NGS_LEN --inner-distance $NGS_INNER --std-dev $NGS_STD --error-rate $NGS_ERR --reads $NGS_READS --fastq1 ${contigs[$i]}.${j}_1.fastq --fastq2 ${contigs[$i]}.${j}_2.fastq
            fi
            # generate 50X TGS
            echo -e "[ Generate TGS data from ${j}th sub population genome of ${contigs[$i]} ]"
            if [ ! -f ${contigs[$i]}.${j}_pacbio.fasta ]; then
                python $prog_path/read_pool-seq_pacbio.py --pg ${contigs[$i]}.$j.fa --tgs-maxl $TGS_MAXL --tgs-minl $TGS_MINL --read-length $TGS_MEANL --tgs-alpha $TGS_ALPHA --tgs-loc $TGS_LOC --tgs-beta $TGS_BETA --error-rate $TGS_ERR --error-fraction $ERR_FRAC --reads $TGS_READS --fasta ${contigs[$i]}.${j}_pacbio.fasta
            fi
            # remove intermediate sub-population genome
            rm ${contigs[$i]}.$j.fa
        done
        cd ..
    fi
done


cat */*_pacbio.fasta > TGS.fasta && samtools faidx TGS.fasta
cat */*_1.fastq > NGS_1.fastq
cat */*_2.fastq > NGS_2.fastq
rm */*fast*


# simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 10000 -R 0 --sub-N 50 --germline-count 500 --avg-somatic-count 20 --min-distance 500 --depth 50
# simulation_protocol.sh -d ./ -r HG02716_template_h1.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 1000 -R 0 --sub-N 20 --germline-count 3000 --avg-somatic-count 20 --min-distance 120000 --depth 50
# simulation_protocol.sh -d ./ -r HG02716_template_h2.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 3000 --avg-somatic-count 20 --min-distance 120000 --depth 50
# simulation_protocol.sh -d ./ -r HG02716_template_h1.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 0 --avg-somatic-count 1000 --min-distance 25000 --depth 5

# for((i=0; i<$contigs_count; i++))
# do
# N_SUB=200
# TGS_READS=`awk -v depth=$DEPTH -v g_l=$genome_size -v c_l=${contigs_size[$i]} -v n_sub=$N_SUB -v r_l=10000 'BEGIN{ratio=c_l/g_l; tgs_reads=int(ratio*(depth*g_l/r_l)/n_sub); print tgs_reads}'`
# echo $TGS_READS
# done

# tmp_ins_count_final=(${INS_NUMs[*]})
# for((i=0; i<$contigs_count; i++))
# do
#     if [ $i -lt $(($contigs_count - 1)) ]; then
#         tmp_ins_count=(`echo ${INS_NUMs[*]} | awk -v r_l=${genome_size} -v q_l=${contigs_size[$i]} 'BEGIN{ratio=q_l/r_l} {for(i=1;i<=NF;i++){print int(ratio*$i+0.5)}}'`)
#         echo $i: ${tmp_ins_count[*]}
        
#         for((j=0; j<${#tmp_ins_count[*]}; j++))
#         do
#             tmp_ins_count_final[$j]=$((${tmp_ins_count_final[$j]} - ${tmp_ins_count[$j]}))
#         done        
#     elif [ $i -eq $(($contigs_count - 1)) ]; then
#         tmp_ins_count=(${tmp_ins_count_final[*]})
#         echo $i: ${tmp_ins_count[*]}
#     fi
# done | cut -f 1

# awk '$1~/^>/{gsub(/-/,"_"); print $0} $1!~/^>/{gsub(/[WSRYKMBDHVN]/,"T"); print $0}' dm3.transposon.fa > dm3.transposon_for_simulaTE.fa && samtools faidx dm3.transposon_for_simulaTE.fa
# nohup simulation_protocol.sh -d ./ -r line_26_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 10000 --sub-N 50 --germline-count 500 --avg-somatic-count 20 --min-distance 500 --depth 10 &

# nohup python /data/tusers/zhongrenhu/Software/bin/change_degenerate_bases.py dm3.transposon.fa dm3.transposon_for_simulaTE.fa &
# awk '$1~/^>/{gsub(/-/,"_"); print $0} $1!~/^>/{print $0}' dm3.transposon_for_simulaTE.fa > t.fa && mv t.fa dm3.transposon_for_simulaTE.fa && samtools faidx dm3.transposon_for_simulaTE.fa
# nohup python /data/tusers/zhongrenhu/Software/bin/change_degenerate_bases.py line_28_template.fa test.fa &
# nohup simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 10000 --sub-N 50 --germline-count 500 --avg-somatic-count 20 --min-distance 500 --depth 50 &
