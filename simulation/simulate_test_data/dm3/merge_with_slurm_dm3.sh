#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH -c 36
#SBATCH --array=1-6
#SBATCH --partition=12hours
#SBATCH --output=./logs/Merge/Dm3/merge-log-dm3-%A-%a.out


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
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/benchmark/merge_filelist
[ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/template/line_28/line_28_template.fa.fai`)


###############################
### Get files from FILELIST ###
###############################
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
SAMPLE=${toks[0]}


###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""


###############################
### Copy inputs to temp dir ###
###############################
echo "Copying input files from ${SAMPLE}..."
for contig in ${CONTIGS[*]}
do
    mkdir ${contig}
    cp ${SAMPLE}/${contig}/TGS.fasta ${contig}
    cp ${SAMPLE}/${contig}/NGS_1.fastq ${contig}
    cp ${SAMPLE}/${contig}/NGS_2.fastq ${contig}
done
echo "Done."
echo ""


###############
### Process ###
###############
echo "Merging..."
cat */TGS.fasta >> TGS.fasta
cat */NGS_1.fastq >> NGS_1.fastq
cat */NGS_2.fastq >> NGS_2.fastq
echo "Done." && echo ""


#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp TGS.fasta $SAMPLE
cp NGS_1.fastq $SAMPLE
cp NGS_2.fastq $SAMPLE
echo "Done."
echo ""


################
### Clean up ###
################
echo "Cleaning up..."
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
#echo "/tmp contents:"
#ls -lh /tmp
echo ""
echo "Script complete."
date