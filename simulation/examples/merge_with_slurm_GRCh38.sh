#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=40G
#SBATCH -c 36
#SBATCH --array=1-150%20
#SBATCH --partition=12hours
#SBATCH --output=./logs/merge-log-%A-%a.out


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
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/v4/merge_filelist
[ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_germ/HG02716_template.fa.fai`)


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
    mkdir ${contig} && cp ${SAMPLE}/${contig}/TGS.fasta ${contig}
done
echo "Done."
echo ""


###############
### Process ###
###############
echo "Merging..."
cat */TGS.fasta >> TGS.fasta
echo "Done." && echo ""


#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp TGS.fasta $SAMPLE
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
