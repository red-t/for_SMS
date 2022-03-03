#! /bin/bash

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-n <name>\tsample name (id)."
    echo -e "\t-s <F/M>\tsample sex, F/M."
    echo -e "\t-f <fasta>\tindexed female reference genome FASTA."
    echo -e "\t-m <fasta>\tindexed male reference genome FASTA."
    echo -e "\t-v <path>\tdirectory to indexed phased SNVs VCFs (1000 genomes for human)."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":n:s:f:m:v:h" OPTION; do
    case $OPTION in
        n)  NAME=$OPTARG;;
        s)  SEX=$OPTARG;;
        f)  F_FASTA=$OPTARG;;
        m)  M_FASTA=$OPTARG;;
        v)  VCF_PATH=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

export PATH=$PATH:/data/tusers/zhongrenhu/Software/anaconda3/bin/
[ ! -z ${F_FASTA} ] && F_CHROMS=($(cut -f 1 ${F_FASTA}.fai))
[ ! -z ${M_FASTA} ] && M_CHROMS=($(cut -f 1 ${M_FASTA}.fai))


### checking ###
echo -e "which VISOR:\t`which VISOR`"
echo -e "NAME:\t${NAME}"
echo -e "SEX:\t${SEX}"
echo -e "VCF_PATH:\t${VCF_PATH}"
echo -e "F_FASTA:\t${F_FASTA}"
echo -e "M_FASTA:\t${M_FASTA}"
echo -e "F_CHROMS:\t${F_CHROMS[*]}"
echo -e "M_CHROMS:\t${M_CHROMS[*]}"
### checking ###


if [ ! -f ${NAME}.tmp.F.${CHROM}.snp.bcf ];then
    mkdir ${NAME}_templateswithsnp

    if [ ${SEX}=="F" ];then
        for CHROM in ${F_CHROMS[*]}
        do
            VCF=`readlink -f ${VCF_PATH%/}/*${CHROM}.*vcf.gz`
            bcftools view --threads 10 -v snps -O b -o ${NAME}.tmp.F.${CHROM}.snp.bcf -s ${NAME} -m2 -M2 -c1 -C1 ${VCF}
            bcftools index ${NAME}.tmp.F.${CHROM}.snp.bcf
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${NAME}.tmp.F.${CHROM}.snp.bcf | grep "1|0" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' >> ${NAME}.tmp.F.snp.h1.bed
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${NAME}.tmp.F.${CHROM}.snp.bcf | grep "0|1" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' >> ${NAME}.tmp.F.snp.h2.bed
        done

        VISOR HACk -g ${F_FASTA} -b ${NAME}.tmp.F.snp.h1.bed ${NAME}.tmp.F.snp.h2.bed -o ${NAME}_templateswithsnp
        mv ${NAME}.tmp.F* ${NAME}_templateswithsnp
    fi

    if [ ${SEX}=="M" ];then
        mkdir haplo2
        for CHROM in ${F_CHROMS[*]}
        do
            VCF=`readlink -f ${VCF_PATH%/}/*${CHROM}.*vcf.gz`
            bcftools view --threads 10 -v snps -O b -o ${NAME}.tmp.M.${CHROM}.snp.bcf -s ${NAME} -m2 -M2 -c1 -C1 ${VCF}
            bcftools index ${NAME}.tmp.M.${CHROM}.snp.bcf
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${NAME}.tmp.M.${CHROM}.snp.bcf | grep "1|0" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' >> ${NAME}.tmp.M.snp.h1.bed
            if [ ${CHROM}!="chrX" ];then
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${NAME}.tmp.M.${CHROM}.snp.bcf | grep "0|1" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' >> ${NAME}.tmp.M.snp.h2.bed
            fi
        done

        VISOR HACk -g ${F_FASTA} -b ${NAME}.tmp.M.snp.h1.bed -o ${NAME}_templateswithsnp
        VISOR HACk -g ${M_FASTA} -b ${NAME}.tmp.M.snp.h2.bed -o haplo2
        mv haplo2/h1.fa haplo2/h2.fa && mv haplo2/h1.fa.fai haplo2/h2.fa.fai && mv haplo2/h2.fa* ${NAME}_templateswithsnp && rm -r haplo2
        mv ${NAME}.tmp.M* ${NAME}_templateswithsnp
    fi
fi  


# HG02716 Female
# HG02610 Male
# /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_X.fa
# /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_Y.fa
# /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GGVP/ALL_GGVP.chrX.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz