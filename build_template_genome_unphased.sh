#! /bin/bash

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-n <name(s)>\tinbred line name(s) (id)."
    echo -e "\t-r <fasta>\tindexed reference genome FASTA."
    echo -e "\t-v <path>\tdirectory to indexed SNVs VCFs with genotype (DGRP for D.mel)."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":n:r:v:h" OPTION; do
    case $OPTION in
        n)  NAMEs=($OPTARG);;
        r)  FASTA=$OPTARG;;
        v)  VCF_PATH=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

export PATH=$PATH:/data/tusers/zhongrenhu/Software/anaconda3/bin/


### checking ###
echo -e "which VISOR:\t`which VISOR`"
echo -e "NAME:\t${NAMEs[*]}"
echo -e "VCF_PATH:\t${VCF_PATH}"
echo -e "FASTA:\t${FASTA}"
### checking ###


for NAME in ${NAMEs[*]}
do
    echo -e "[ BUILD TEMPLATE GENOME FROM ${NAME} START. ]"
    if [ ! -f ${NAME}.tmp.snp.bed ];then
        [ ! -d ${NAME}_templateswithsnp ] && mkdir ${NAME}_templateswithsnp
        VCF=`readlink -f ${VCF_PATH%/}/*vcf.gz`
        bcftools view --threads 10 -v snps -O b -o ${NAME}.tmp.snp.bcf -s ${NAME} -m2 -M2 -c1 -C2 ${VCF}
        bcftools index ${NAME}.tmp.snp.bcf
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${NAME}.tmp.snp.bcf | grep '1/1' | awk 'OFS=FS="\t"''{print "chr"$1, ($2 -1), $2, "SNP", $4, "0"}' >> ${NAME}.tmp.snp.bed

        VISOR HACk -g ${FASTA} -b ${NAME}.tmp.snp.bed -o ${NAME}_templateswithsnp
        mv ${NAME}.tmp.* ${NAME}_templateswithsnp && mv ${NAME}_templateswithsnp/h1.fa ${NAME}_templateswithsnp/${NAME}_template.fa && mv ${NAME}_templateswithsnp/h1.fa.fai ${NAME}_templateswithsnp/${NAME}_template.fa.fai
    fi  

    if [ -f ${NAME}.tmp.snp.bed ];then
        VISOR HACk -g ${FASTA} -b ${NAME}.tmp.snp.bed -o ${NAME}_templateswithsnp
        mv ${NAME}.tmp.* ${NAME}_templateswithsnp && mv ${NAME}_templateswithsnp/h1.fa ${NAME}_templateswithsnp/${NAME}_template.fa && mv ${NAME}_templateswithsnp/h1.fa.fai ${NAME}_templateswithsnp/${NAME}_template.fa.fai
    fi

    echo -e "[ BUILD TEMPLATE GENOME FROM ${NAME} FINISH. ]"
done



# /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.fa
# /data/tusers.ds/zhongrenhu/for_SMS/reference/dm3/DGRP/

# line_21 line_26 line_28 line_31 line_32 line_38 line_40 line_41 line_42 line_45
# line_48 line_49 line_57 line_59 line_69 line_73 line_75 line_83 line_85 line_88
# line_91 line_93 line_100 line_101 line_105 line_109 line_129 line_136 line_138 line_142
# line_149 line_153 line_158 line_161 line_176 line_177 line_181 line_189 line_195 line_208
# line_217 line_223 line_227 line_228 line_229 line_233 line_235 line_237 line_239 line_256
# line_280 line_287 line_301 line_303 line_304 line_306 line_307 line_309 line_310 line_313
# line_315 line_317 line_318 line_319 line_320 line_321 line_324 line_325 line_332 line_335
# line_336 line_338 line_340 line_348 line_350 line_352 line_354 line_355 line_356 line_357
# line_358 line_359 line_360 line_361 line_362 line_365 line_367 line_370 line_371 line_373
# line_374 line_375 line_377 line_379 line_380 line_381 line_382 line_383 line_385 line_386
# line_390 line_391 line_392 line_395 line_397 line_399 line_405 line_406 line_409 line_426