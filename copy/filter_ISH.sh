#!/bin/bash

# filterring reads (alignments) with insert fragment(s) longer than specified length.
# insert fragments: fragments correspond to I, S, H CIGAR OPERATION.


######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-b <bam>\talignment in BAM format, with samtools index first."
    echo -e "\t-@ <cores>\ttCPU cores numbers."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":b:@:l:h" OPTION; do
    case $OPTION in
        b)  BAM=$OPTARG;;
        @)  CORE=$OPTARG;;
        l)  LENGTH=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


######## Processing ########
$PREFIX=`basename ${BAM%.bam}`
echo -e "filterring reads (alignments) with insert fragment(s) longer ${LENGTH} for ${PREFIX}.bam: "
start=$(date +%s)

# Get BAM file header
samtools view -H ${BAM} > tmp.header

# filterring reads parallelly
samtools view -@ ${CORE} ${BAM} > tmp.sam
parallel --pipepart -a tmp.sam --block 50M -j ${CORE} -q awk -v len=${LENGTH} '{split($6, a, /[ISH]/); max=0; for(i in a){match(a[i], /([0-9]+$)/, b); if(b[1]>max){max=b[1]}}; if(max>=len){print $0}}' | cat tmp.header - | samtools view -@ ${CORE} -bhS - | samtools sort -@ ${CORE} -o ${PREFIX}_filtered.bam - && samtools index -@ ${CORE} ${PREFIX}_filtered.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds.