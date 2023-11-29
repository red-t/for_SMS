#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -c 20
#SBATCH --array=1-3450%20
#SBATCH --partition=12hours
#SBATCH --output=./logs/simulation-log-%A-%a.out


#########################################
### Print a little info for debugging ###
#########################################
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


################################
### Parameter initialization ###
################################
# [ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/v1/HG02716_0/simulation_filelist
# [ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/v3/HG02716_0/simulation_filelist
# [ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/v4_test/HG02716_0_ccs/simulation_filelist
[ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/v4/HG02716_0_ccs/simulation_filelist

[ -z $NGS_LEN ] && NGS_LEN=150
[ -z $NGS_INNER ] && NGS_INNER=200
[ -z $NGS_STD ] && NGS_STD=20
[ -z $NGS_ERR ] && NGS_ERR=0.0005
[ -z $TGS_MAXL ] && TGS_MAXL=300000
[ -z $TGS_MINL ] && TGS_MINL=100

### Get parameters from FILELIST ###
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
WORKDIR=${toks[0]} && CONTIG=`basename ${WORKDIR}`
TE_FA=${toks[1]}
GENOME_SIZE=${toks[2]}
CONTIG_SIZE=${toks[3]}
POP_SIZE=${toks[4]}
SUB_POP_SIZE=${toks[5]}
DEPTH=${toks[6]}
PROTOCOL=${toks[7]}
[ $PROTOCOL == "ccs" ] && TGS_MEANL=13490
[ $PROTOCOL == "clr" ] && TGS_MEANL=7896
[ $PROTOCOL == "ont" ] && TGS_MEANL=7170
PROG_PATH="/data/tusers/zhongrenhu/for_SMS/bin/simulation"


###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)

### Copy inputs to temp dir ###
echo "Copying ${WORKDIR} to temp directory..."
cp -r $WORKDIR .
echo "Done."
echo ""


##################################
### Activate conda environment ###
##################################
source /zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh
conda activate simulation


########################
### Process the data ###
########################
cd $CONTIG
N_SUB=`awk -v sub_pop=$SUB_POP_SIZE -v pop=$POP_SIZE 'BEGIN{n_sub=int(pop/sub_pop); print n_sub}'`
NGS_READS=`awk -v depth=$DEPTH -v g_l=$GENOME_SIZE -v c_l=$CONTIG_SIZE -v n_sub=$N_SUB -v r_l=150 -v inner=$NGS_INNER 'BEGIN{ratio=c_l/g_l; ngs_reads=int(ratio*(depth*g_l/(2*r_l+inner))/n_sub); print ngs_reads}'`
TGS_READS=`awk -v depth=$DEPTH -v g_l=$GENOME_SIZE -v c_l=$CONTIG_SIZE -v n_sub=$N_SUB -v r_l=$TGS_MEANL 'BEGIN{ratio=c_l/g_l; tgs_reads=int(ratio*(depth*g_l/r_l)/n_sub); print tgs_reads}'`

if [ -f $CONTIG.0.pgd ]; then
    for((j=0; j<$N_SUB; j++))
    do
        # build population geneome
        echo -e "[ Build ${j}th sub population genome of $CONTIG ]"
        if [ ! -f $CONTIG.$j.fa ]; then
            # python $PROG_PATH/build-population-genome.py --pgd $CONTIG.$j.pgd --te-seqs $TE_FA --chassis $CONTIG.tmp.chasis.fasta --output $CONTIG.$j.fa --ins-seq $CONTIG.ins.sequence --sub_idx $j --sub_size $SUB_POP_SIZE
            python $PROG_PATH/build-population-genome.py --pgd $CONTIG.$j.pgd --te-seqs $TE_FA --chassis $CONTIG.tmp.chasis.fasta --output $CONTIG.$j.fa --ins-seq "" --sub_idx $j --sub_size $SUB_POP_SIZE
        fi

        # generate NGS
        # echo -e "[ Generate NGS data from ${j}th sub population genome of $CONTIG ]"
        # if [ ! -f $CONTIG.${j}_2.fastq ]; then
        #    python $PROG_PATH/generate_NGS.py --pg $CONTIG.$j.fa --read-length $NGS_LEN --inner-distance $NGS_INNER --std-dev $NGS_STD --error-rate $NGS_ERR --reads $NGS_READS --fastq1 $CONTIG.${j}_1.fastq --fastq2 $CONTIG.${j}_2.fastq
        # fi

        # generate TGS
        echo -e "[ Generate TGS data from ${j}th sub population genome of $CONTIG ]"
        if [ ! -f $CONTIG.${j}_tgs.fasta ]; then
            python $PROG_PATH/generate_TGS.py --pg $CONTIG.$j.fa --reads $TGS_READS --fasta $CONTIG.${j}_tgs.fasta --tgs-maxl $TGS_MAXL --tgs-minl $TGS_MINL --protocol $PROTOCOL
        fi

        # remove intermediate sub-population genome
        rm $CONTIG.$j.fa
    done
fi

cd ..
[ -f $CONTIG/$CONTIG.0_tgs.fasta ] && cat */*_tgs.fasta > TGS.fasta
# [ -f $CONTIG/$CONTIG.0_1.fastq ] && cat */*_1.fastq > NGS_1.fastq
# [ -f $CONTIG/$CONTIG.0_2.fastq ] && cat */*_2.fastq > NGS_2.fastq


################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp -r TGS* $WORKDIR
# cp -r NGS* $WORKDIR
echo "Done."
echo ""


################
### Clean up ###
################
echo "Cleaning up..."
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
echo ""
echo "Script complete."
date