#!/bin/bash


help_info(){
	
echo -e "\033[40;32;1m
$0
	-g genome_fa !important
    -G genome minimap2 index

    -e transposon fasta file with subfamily  !important
    -E transposon minimap2 index

    -r repeatmasker annotation file

	-i data path of fasta file
    -I data path of bam file

	-o out path to save your reasult

    -f flanksize

	-p prefix
	-c cpu number
	\033[0m"
}

if [ $# -lt 1 ];then
	help_info && exit 1
fi

while getopts "hg:G:e:E:r:i:I:o:f:p:c:" OPTION;
do
	case $OPTION in
		h)	help_info && exit 1;;
		g)	GENOME_FA=${OPTARG};;
        G)  GENOME_INDEX=${OPTARG};;
        e)  TRANSPOSON_FA=${OPTARG};;
        E)  TRANSPOSON_INDEX=${OPTARG};;
		i)  FASTA_FILE=${OPTARG};;
        I)  BAM_FILE=${OPTARG};;
        r)  REPEATMASKER_FILE=${OPTARG};;
		o)	OUT_PATH=${OPTARG};;
        f)  FLANKSIZE=${OPTARG};;
        p)  PREFIX=${OPTARG};;
		c)	CPU=${OPTARG};;
	esac
done

# out path
if [ -z ${OUT_PATH} ];then
    OUT_PATH="./SMS"
fi
[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
[ ! -d ${OUT_PATH}/temp ] && mkdir ${OUT_PATH}/temp
[ ! -d mkdir ${OUT_PATH}/anno ] && mkdir mkdir ${OUT_PATH}/anno

# reference file ; genome file and transposon
if [ -z ${GENOME_FA} ];then
    echo "Input Option Error: Please provide the reference genome file"
    help_info && exit 1
elif [ ! -f ${GENOME_FA} ];then
    echo "Input Option Error: There is no ${GENOME_FA}"
    help_info && exit 1
fi

if [ -z ${TRANSPOSON_FA} ];then
    echo "Input Option Error: Please provide the transposon consensus sequence file"
    help_info && exit 1
elif [ ! -f ${TRANSPOSON_FA} ];then
    echo "Input Option Error: There is no ${TRANSPOSON_FA}"
    help_info && exit 1
else
    samtools faidx ${TRANSPOSON_FA} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2}' > ${OUT_PATH}/anno/transposon.size
    TRANSPOSON_SIZE=${OUT_PATH}/anno/transposon.size
fi



# idnex; genome and tansposon index

if [ -z ${GENOME_INDEX} ];then
    echo "Running: minimap2 -d ${OUT_PATH}/anno/genome.mmi ${GENOME_FA}"
    minimap2 -d ${OUT_PATH}/anno/genome.mmi ${GENOME_FA}
    GENOME_INDEX=${OUT_PATH}/anno/genome.mmi
elif [ ! -f ${GENOME_INDEX} ];then
    echo "Input Option Error: There is no ${GENOME_INDEX}"
    help_info && exit 1
fi


if [ -z ${TRANSPOSON_INDEX} ];then
    echo " Running: 'minimap2 -d ${OUT_PATH}/anno/transposon.mmi ${TRANSPOSON_FA} ' "
    minimap2 -d ${OUT_PATH}/anno/transposon.mmi ${TRANSPOSON_FA}
    TRANSPOSON_INDEX=${OUT_PATH}/anno/transposon.mmi
elif [ ! -f ${TRANSPOSON_INDEX} ];then
    echo "Input Option Error: There is no ${TRANSPOSON_INDEX}"
    help_info && exit 1
fi

# data
if [ -z ${FASTA_FILE} ] && [ -z ${BAM_FILE} ];then
    echo "Input Option Error: Please provide the reads fasta file or bam file"
    help_info && exit 1
elif [ ! -f ${FASTA_FILE} ];then
    echo "Input Option Error: There is no ${FASTA_FILE}"
    ehelp_info && exit 1
elif [ ! -f ${BAM_FILE} ];then
    echo "Input Option Error: There is no ${BAM_FILE}"
    help_info && exit 1
fi


if [ -z ${REPEATMASKER_FILE} ];then
    echo "Input Option Error: Please provide repeatmasker file"
    help_info && exit 1
elif [ ! -f ${REPEATMASKER_FILE} ];then
    echo "Input Option Error: There is no ${REPEATMASKER_FILE} "
    help_info && exit 1
fi



# flanksize
if [ -z ${FLANKSIZE} ];then
    echo "Set flanksize = 250 "
    FLANKSIZE=250
fi

# prefix
if [ -z ${PREFIX} ];then
    if [ -z ${FASTA_FILE} ];then
        PREFIX_temp=${BAM_FILE##*/}
        PREFIX=${PREFIX_temp%.*}
    elif [ -z ${BAM_FILE} ];then
        PREFIX_temp=${FASTA_FILE##*/}
        PREFIX=${PREFIX_temp%.*}
    fi
fi
echo ${PREFIX}


# CPU
if [ -z ${CPU} ];then
    echo "Set Cpu num to 10"
    CPU=10
fi



if [ -z ${BAM_FILE} ];then
    minimap2 -t ${CPU} -aYx map-ont ${GENOME_INDEX} ${FASTA_FILE} > ${OUT_PATH}/${PREFIX}.genome.sam
    samtools sort -@ ${CPU} -O bam -o ${OUT_PATH}/${SAMPLE}.genome.bam ${OUT_PATH}/${SAMPLE}.genome.sam
    samtools index ${OUT_PATH}/${SAMPLE}.genome.bam
    BAM_FILE=${OUT_PATH}/${SAMPLE}.genome.bam
fi



if [ ! -f ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ];then
    samtools view -bhSq 0 ${OUT_PATH}/${SAMPLE}.genome.bam > ${OUT_PATH}/${SAMPLE}.q0.bam
    samtools sort -@ ${CPU} -o ${OUT_PATH}/${SAMPLE}.q0.sorted.bam ${OUT_PATH}/${SAMPLE}.q0.bam
    samtools index -@ ${CPU} ${OUT_PATH}/${SAMPLE}.q0.sorted.bam
    BAM_FILE=${OUT_PATH}/${SAMPLE}.q0.sorted.bam
fi


# python lrft2.py ${OUT_PATH}/${SAMPLE}_pacbio.${dep}X.q0.region.sorted.bam ${OUT_PATH}/SMS ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.mmi ${ANNO_PATH}/dm3/dm3.transposon_for_simulaTE.size 300

python lrft2.py ${BAM_FILE} ${OUT_PATH} ${TRANSPOSON_INDEX} ${TRANSPOSON_SIZE} ${FLANKSIZE} ${PREFIX} ${GENOME_FA} ${GENOME_INDEX} ${REPEATMASKER_FILE}




