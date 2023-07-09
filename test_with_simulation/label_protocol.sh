#! /bin/bash

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-b <bam>"
    echo -e "\t-d <workdir>"
    echo -e "\t-r <repeat.bed>"
    echo -e "\t-g <gap.bed>"
    echo -e "\t-T <ref_te.fa>"
    echo -e "\t-x <preset>"
    echo -e "\t-p <nprocess>"
    echo -e "\t-t <nthreads>"
    echo -e "\t-l <minl>"
    echo -e "\t-L <maxdist>"
    echo -e "\t-G <path> path of TEMP3.py"
    echo -e "\t-s <subsize> size of sub-population genome"
    echo -e "\t-i <reftid> specify a tid, which will be used to estimate background information"
    echo -e "\t-h \tShow this information"
}


######## Getting Parameters ########
ARGS=`getopt -o b:d:r:g:T:x:p:t:l:L:G:s:i:h -n "$0" -- "$@"`
if [ $? != 0 ]; then
    echo "Terminating..."
    exit 1
fi

eval set -- "${ARGS}"

while true
do
    case "$1" in
        -b) 
            BAM=$2
            shift 2
            ;;
        -d) 
            WORK_DIR=$2
            shift 2
            ;;
        -r)
            REPEAT=$2
            shift 2
            ;;
        -g)
            GAP=$2
            shift 2
            ;;
        -T)
            TE=$2
            shift 2
            ;;
        -x)
            PRESET=$2
            shift 2
            ;;
        -p)
            NPROCESS=$2
            shift 2
            ;;
        -t)
            NTHREADS=$2
            shift 2
            ;;
        -l)
            MINL=$2
            shift 2
            ;;
        -L)
            MAXDIST=$2
            shift 2
            ;;
        -G)
            TEMP3_PATH=$2
            shift 2
            ;;
        -s)
            SUBSIZE=$2
            shift 2
            ;;
        -i)
            REFTID=$2
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



[ -z $SUBSIZE ] && SUBSIZE=50
[ -z $REFTID ] && REFTID=0
### Checking Parameters ###
echo -e "label_protocol.sh parameters:"
echo -e "BAM:\t${BAM}"
echo -e "WORK_DIR:\t${WORK_DIR}"
echo -e "REPEAT:\t${REPEAT}"
echo -e "GAP:\t${GAP}"
echo -e "TE:\t${TE}"
echo -e "PRESET:\t${PRESET}"
echo -e "NPROCESS:\t${NPROCESS}"
echo -e "NTHREADS:\t${NTHREADS}"
echo -e "MINL:\t${MINL}"
echo -e "MAXDIST:\t${MAXDIST}"
echo -e "TEMP3_PATH:\t${TEMP3_PATH}"
echo -e "SUBSIZE:\t${SUBSIZE}"


### Generate All Candidates ###
PROG_PATH=`readlink -f $0` && PROG_PATH=`dirname ${PROG_PATH}`
python ${PROG_PATH}/../develop/TEMP3.py -b ${BAM} -r ${REPEAT} -g ${GAP} -T ${TE} -p ${NPROCESS} -t ${NTHREADS} -l ${MINL} -i ${REFTID}

# merge candidate clusters & segments
cat tmp_clt_*txt > tmp_all_clt.txt && cut -f 1-6 tmp_all_clt.txt > tmp_all_candidates.bed && rm tmp_clt_*txt
cat tmp_seg_*txt > tmp_all_seg.txt && rm tmp_seg_*txt

# merge candidate alignments
samtools view -H tmp_candidates_alignments.0.bam > header
for bam in tmp_candidates_alignments*bam
do
    samtools view -@ ${NTHREADS} ${bam} >> tmp.sam
done
cat header tmp.sam | samtools view -@ ${NTHREADS} -bhS - | samtools sort -@ ${NTHREADS} -o tmp_all_candidates_alignments.bam - && samtools index tmp_all_candidates_alignments.bam
rm tmp.sam && rm tmp_candidates_alignments*bam && rm header


### Filtering ###
cut -f 1-6 ${WORK_DIR}/*/*ins.summary > AllIns.bed && bgzip AllIns.bed && tabix AllIns.bed.gz
python ${PROG_PATH}/Filter_TP_and_FP.py --workdir ${WORK_DIR} --bam tmp_all_candidates_alignments.bam --qbed tmp_all_candidates.bed --rbed AllIns.bed.gz --segf tmp_all_seg.txt --subsize ${SUBSIZE}
python ${PROG_PATH}/rename.py
${PROG_PATH}/merge_TP_and_FP.sh tmp_all_candidates_alignments.bam


### Clean Up###
rm fp_*bam* && rm tp_*bam*
rm tmp*fa* && rm tmp*bam*
rm tmp_all_candidates.bed
rm tmp_all_clt.txt && rm tmp_all_seg.txt
rm header
