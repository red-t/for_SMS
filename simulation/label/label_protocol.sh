#! /bin/bash

### Help Information ###
function help_info(){
    echo `basename $0`
    echo -e "\t-b <bam>"
    echo -e "\t-d <workDir> path to simulation workdir"
    echo -e "\t-r <repeat.bed>"
    echo -e "\t-g <gap.bed>"
    echo -e "\t-P <numProcess>"
    echo -e "\t-T <numThread>"
    echo -e "\t--refTe <refTe.fa>"
    echo -e "\t--blacklist <blacklist.fa>"
    echo -e "\t--minLen <minLen>"
    echo -e "\t--TEMP3 <path> path to TEMP3.py"
    echo -e "\t--subSize <subSize> size of sub-population genome"
    echo -e "\t--germModel <path> path to germline insertion model"
    echo -e "\t--somaModel <path> path to somatic insertion model"
    echo -e "\t-h \tShow this information"
}


### Getting Parameters ###
ARGS=`getopt -o b:d:r:g:P:T:h --long refTe:,blacklist:,minLen:,TEMP3:,subSize:,germModel:,somaModel: -n "$0" -- "$@"`
if [ $? != 0 ]; then
    echo "Terminating..."
    exit 1
fi

eval set -- "${ARGS}"

while true
do
    case "$1" in
        -b) 
            bam=$2
            shift 2
            ;;
        -d) 
            workDir=$2
            shift 2
            ;;
        -r)
            repeatBed=$2
            shift 2
            ;;
        -g)
            gapBed=$2
            shift 2
            ;;
        -P)
            numProcess=$2
            shift 2
            ;;
        -T)
            numThread=$2
            shift 2
            ;;
        --refTe)
            refTe=$2
            shift 2
            ;;
        --blacklist)
            blackList=$2
            shift 2
            ;;
        --minLen)
            minLen=$2
            shift 2
            ;;
        --TEMP3)
            TEMP3=$2
            shift 2
            ;;
        --subSize)
            subSize=$2
            shift 2
            ;;
        --germModel)
            germModel=$2
            shift 2
            ;;
        --somaModel)
            somaModel=$2
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


### Check Parameters ###
echo -e "label_protocol.sh parameters:"
echo -e "BAM:\t${bam}"
echo -e "WORK_DIR:\t${workDir}"
echo -e "REPEAT:\t${repeatBed}"
echo -e "GAP:\t${gapBed}"
echo -e "TE:\t${refTe}"
echo -e "BLACKLIST:\t${blacklist}"
echo -e "NPROCESS:\t${numProcess}"
echo -e "NTHREADS:\t${numThread}"
echo -e "MINL:\t${minLen}"
echo -e "TEMP3_PATH:\t${TEMP3}"
echo -e "SUBSIZE:\t${subSize}"


### Labeling on blades ###
programPath=`readlink -f $0` && programPath=`dirname ${programPath}`
if [ ! -f 'tmp_all_seg.txt' ]; then
    ### Generate All Candidates ###
    python ${TEMP3}/TEMP3.py \
        --bam ${bam} \
        --repeat ${repeatBed} \
        --gap ${gapBed} \
        --blacklist ${blackList} \
        --refte ${refTe} \
        --germ ${germModel} \
        --soma ${somaModel} \
        --numprocess ${numProcess} \
        --numthread ${numThread} \
        --minSegLen ${minLen}
    
    ### Merge candidate clusters & segments ###
    cat tmp_clt_*txt > tmp_all_clt.txt && cut -f 1-6 tmp_all_clt.txt > tmp_all_candidates.bed && rm tmp_clt_*txt
    cat tmp_seg_*txt > tmp_all_seg.txt && rm tmp_seg_*txt
    
    ### Merge candidate alignments ###
    samtools view -H tmp_candidates_alignments.0.bam > header
    for tmpBam in tmp_candidates_alignments*bam
    do
        samtools view -@ ${numThread} ${tmpBam} >> tmp.sam
    done

    cat header tmp.sam \
        | samtools view -@ ${numThread} -bhS - \
        | samtools sort -@ ${numThread} -o tmp_all_candidates_alignments.bam - \
        && samtools index tmp_all_candidates_alignments.bam
    
    rm tmp.sam && rm tmp_candidates_alignments*bam && rm header
else
    ### Filtering ###
    cut -f 1-6 ${workDir}/*/*ins.summary > AllIns.bed && bgzip AllIns.bed && tabix AllIns.bed.gz
    python ${programPath}/Filter_TP_and_FP.py \
        --workdir ${workDir} \
        --bam tmp_all_candidates_alignments.bam \
        --qbed tmp_all_candidates.bed \
        --rbed AllIns.bed.gz \
        --segf tmp_all_seg.txt \
        --subsize ${subSize}
    python ${programPath}/rename.py

    ${programPath}/merge_TP_and_FP.sh tmp_all_candidates_alignments.bam

    ### Clean Up###
    rm fp_*bam* && rm tp_*bam*
    rm tmp*fa* && rm tmp*bam*
    rm tmp_all_candidates.bed
    rm tmp_all_clt.txt && rm tmp_all_seg.txt
    rm header
fi
