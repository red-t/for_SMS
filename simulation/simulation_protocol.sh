#! /bin/bash

########################
### Help Information ###
########################
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
    echo -e "\t--species <human/fly>\tthe transposon library to the corresponding species will be used"
    echo -e "\t--protocol <ccs/clr/ont>\tthe error rate, error fraction and length distribution to the corresponding protocol will be used"
    echo -e "\t--truncProb <float>\tProbability of constructing a truncation on an inserted sequence"
    echo -e "\t--nestProb <float>\tProbability of constructing a nested insertion on an inserted sequence"
    echo -e "\t--mode <1/2>\tGenerate data for: (1)model training or (2)benchmarking"
    echo -e "\t-h \tShow this information"
}


##########################
### Getting Parameters ###
##########################
ARGS=`getopt -o d:r:t:N:R:h --long sub-N:,germline-count:,avg-somatic-count:,min-distance:,depth:,ngs-len:,ngs-inner:,ngs-std:,ngs-err:,species:,protocol:,truncProb:,nestProb:,mode: -n "$0" -- "$@"`
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
        --species)
            SPECIES=$2
            shift 2
            ;;
        --protocol)
            PROTOCOL=$2
            shift 2
            ;;
        --truncProb)
            truncProb=$2
            shift 2
            ;;
        --nestProb)
            nestProb=$2
            shift 2
            ;;
        --mode)
            mode=$2
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


##########################
### Default Parameters ###
##########################
[ -z $D_RATE ] && D_RATE=0
[ -z $NGS_LEN ] && NGS_LEN=150
[ -z $NGS_INNER ] && NGS_INNER=200
[ -z $NGS_STD ] && NGS_STD=20
[ -z $NGS_ERR] && NGS_ERR=0.0005

[ -z $SPECIES ] && SPECIES="fly"
[ -z $PROTOCOL ] && PROTOCOL="ccs"
[ -z $TGS_MINL ] && TGS_MINL=100
[ -z $TGS_MAXL ] && TGS_MAXL=300000
[ $PROTOCOL == "ccs" ] && TGS_MEANL=13490
[ $PROTOCOL == "clr" ] && TGS_MEANL=7896
[ $PROTOCOL == "ont" ] && TGS_MEANL=7170


###########################
### Checking Parameters ###
###########################
echo -e "Parameters:"
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
echo -e "SPECIES:\t${SPECIES}"
echo -e "PROTOCOL:\t${PROTOCOL}"
echo -e "TGS_MINL:\t${TGS_MINL}"
echo -e "TGS_MAXL:\t${TGS_MAXL}"
echo -e "TGS_MEANL:\t${TGS_MEANL}"


####################
### Program Path ###
####################
source /zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh
conda activate simulation
prog_path=`readlink -f $0` && prog_path=`dirname $prog_path`


#############################################
### Generate Population Genome Definition ###
#############################################
S_COUNT=$(($AVG_S_COUNT * $POP_SIZE))
genome_size=`cut -f 2 ${REF_FA}.fai | awk 'BEGIN{n=0} {n=n+$1} END{print n}'`
contigs=(`cut -f 1 ${REF_FA}.fai`)
contigs_size=(`cut -f 2 ${REF_FA}.fai`)
contigs_count=`wc -l ${REF_FA}.fai | awk '{print $1}'`
delta_g=0
delta_s=0
for((i=0; i<$contigs_count; i++))
do
    # Calculate insertions numbers of each contig based on their length
    if [ $i -lt $(($contigs_count - 1)) ]; then
        tmp_g_count=`awk -v g_count=${G_COUNT} -v g_l=${genome_size} -v c_l=${contigs_size[$i]} 'BEGIN{ratio=c_l/g_l; print int(g_count*ratio+0.5)}'`
        tmp_s_count=`awk -v s_count=${S_COUNT} -v g_l=${genome_size} -v c_l=${contigs_size[$i]} 'BEGIN{ratio=c_l/g_l; print int(s_count*ratio+0.5)}'`
        delta_g=$(($delta_g + $tmp_g_count))
        delta_s=$(($delta_s + $tmp_s_count))
    elif [ $i -eq $(($contigs_count - 1)) ]; then
        tmp_g_count=$(($G_COUNT - $delta_g))
        tmp_s_count=$(($S_COUNT - $delta_s))
    fi
    
    # Extract contig sequence AND generate header & body of the pgd-file
    echo -e "Defining TE landscapes for ${contigs[$i]}:"
    if [ ! -f ${contigs[$i]}/${contigs[$i]}.ins.summary ]; then
        samtools faidx ${REF_FA} ${contigs[$i]} > ${contigs[$i]}.tmp.chasis.fasta && samtools faidx ${contigs[$i]}.tmp.chasis.fasta
        python $prog_path/define_population_genome.py --chassis ${contigs[$i]}.tmp.chasis.fasta --te-seqs ${TE_FA} --N ${POP_SIZE} --divergence-rate ${D_RATE} --germline-count ${tmp_g_count} --somatic-count ${tmp_s_count} --min-distance ${MIN_DIST} --species ${SPECIES} --truncProb ${truncProb} --nestProb ${nestProb} --mode ${mode}
        
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

        mkdir ${contigs[$i]} && mv ${contigs[$i]}.* ${contigs[$i]}
    else
        echo -e "The population-genome-definition (pgd) files of ${contigs[$i]} already exist, skip."
    fi
done


##########################
### Generate Raw Reads ###
##########################
# for((i=0; i<$contigs_count; i++))
# do
#     if [ -f ${contigs[$i]}/${contigs[$i]}.0.pgd ]; then
#         cd ${contigs[$i]}
#         N_SUB=`awk -v sub_pop=$SUB_POP_SIZE -v pop=$POP_SIZE 'BEGIN{n_sub=int(pop/sub_pop); print n_sub}'`
#         NGS_READS=`awk -v depth=$DEPTH -v g_l=$genome_size -v c_l=${contigs_size[$i]} -v n_sub=$N_SUB -v r_l=150 'BEGIN{ratio=c_l/g_l; ngs_reads=int(ratio*(depth*g_l/(2*r_l))/n_sub); print ngs_reads}'`
#         TGS_READS=`awk -v depth=$DEPTH -v g_l=$genome_size -v c_l=${contigs_size[$i]} -v n_sub=$N_SUB -v r_l=$TGS_MEANL 'BEGIN{ratio=c_l/g_l; tgs_reads=int(ratio*(depth*g_l/r_l)/n_sub); print tgs_reads}'`
#         for((j=0; j<$N_SUB; j++))
#         do
#             # build population geneome
#             echo -e "[ Build ${j}th sub population genome of ${contigs[$i]} ]"
#             if [ ! -f ${contigs[$i]}.$j.fa ]; then
#                 # python $prog_path/simulate/build-population-genome.py --pgd ${contigs[$i]}.$j.pgd --te-seqs $TE_FA --chassis ${contigs[$i]}.tmp.chasis.fasta --output ${contigs[$i]}.$j.fa --ins-seq ${contigs[$i]}.ins.sequence --sub_idx $j --sub_size $SUB_POP_SIZE
#                 python $prog_path/build-population-genome.py --pgd ${contigs[$i]}.$j.pgd --te-seqs $TE_FA --chassis ${contigs[$i]}.tmp.chasis.fasta --output ${contigs[$i]}.$j.fa --ins-seq "" --sub_idx $j --sub_size $SUB_POP_SIZE
#             fi
#             # generate NGS
#             echo -e "[ Generate NGS data from ${j}th sub population genome of ${contigs[$i]} ]"
#             if [ ! -f ${contigs[$i]}.${j}_2.fastq ]; then
#                 python $prog_path/generate_NGS.py --pg ${contigs[$i]}.$j.fa --read-length $NGS_LEN --inner-distance $NGS_INNER --std-dev $NGS_STD --error-rate $NGS_ERR --reads $NGS_READS --fastq1 ${contigs[$i]}.${j}_1.fastq --fastq2 ${contigs[$i]}.${j}_2.fastq
#             fi
#             # generate TGS
#             echo -e "[ Generate TGS data from ${j}th sub population genome of ${contigs[$i]} ]"
#             if [ ! -f ${contigs[$i]}.${j}_pacbio.fasta ]; then
#                 python $prog_path/generate_TGS.py --pg ${contigs[$i]}.$j.fa --reads $TGS_READS --fasta ${contigs[$i]}.${j}_pacbio.fasta --tgs-maxl $TGS_MAXL --tgs-minl $TGS_MINL --protocol ${PROTOCOL}
#             fi
#             # remove intermediate sub-population genome
#             rm ${contigs[$i]}.$j.fa
#         done
#         cd ..
#     fi
# done





##############
### dm3 v1 ###
##############
# simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 10000 -R 0 --sub-N 50 --germline-count 500 --avg-somatic-count 20 --min-distance 500 --depth 50

##############
### dm3 v3 ###
##############
# simulation_protocol.sh -d ./ -r line_21_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 100 --min-distance 6000 --depth 50

###################
### dm3 v4 test ###
###################
# simulation_protocol.sh -d ./ -r line_21_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 500 --avg-somatic-count 20 --min-distance 9000 --depth 50 --species fly --protocol ccs

##############
### dm3 v4 ###
##############
# simulation_protocol.sh -d ./ -r line_21_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 50 --min-distance 9000 --depth 50 --species fly --protocol ccs --truncProb 0.1 --nestProb 0.1 --mode 1



#########################
### GRCh38 v1 somatic ###
#########################
# simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 1000 -R 0 --sub-N 20 --germline-count 0 --avg-somatic-count 20 --min-distance 4000 --depth 50

##########################
### GRCh38 v1 germline ###
##########################
# simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 5000 --avg-somatic-count 0 --min-distance 100000 --depth 50

#################
### GRCh38 v3 ###
#################
# simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 5000 --avg-somatic-count 1000 --min-distance 4000 --depth 50

######################
### GRCh38 v4 test ###
######################
# simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 2000 --avg-somatic-count 20 --min-distance 9000 --depth 50 --species human --protocol ccs

#################
### GRCh38 v4 ###
#################
# simulation_protocol.sh -d ./ -r HG02716_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 3000 --avg-somatic-count 500 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.1 --nestProb 0.1 --mode 1