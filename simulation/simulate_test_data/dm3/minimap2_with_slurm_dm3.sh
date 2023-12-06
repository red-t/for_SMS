#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH -c 36
#SBATCH --array=1-6
#SBATCH --partition=12hours
#SBATCH --output=./logs/Minimap2/Dm3/minimap2-log-%A-%a.out


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
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/benchmark/minimap2_filelist
[ -z $SAM_FORMAT ] && SAM_FORMAT=1
[ -z $THREADS ] && THREADS=30


###############################
### Get files from FILELIST ###
###############################
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
REF=${toks[0]}
QUERY=${toks[1]}
PRESET=${toks[2]}


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
echo "Copying input file ${QUERY} to temp directory..."
cp $REF .
cp $QUERY .
echo "Done."
echo ""


####################################################
### Get filenames in temp by cutting of the path ###
####################################################
PREFIX=`basename ${QUERY%.f*a*}` && QUERY_FA=`basename $QUERY` && REF_FA=`basename $REF`
[ -z $OUTPUTDIR ] && OUTPUTDIR=`dirname $QUERY`


##################################
### Activate conda environment ###
##################################
source /zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh
conda activate TEMP3


########################
### Process the data ###
########################
echo "Running minimap2 for ${PREFIX}"
if [ -n $SAM_FORMAT ];then
    minimap2 -aYx $PRESET --MD -t $THREADS $REF_FA $QUERY_FA | samtools view -@ $THREADS -bhS - | samtools sort -@ $THREADS -o $PREFIX.bam -
    samtools index -@ $THREADS $PREFIX.bam
fi
echo ""


#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
mkdir ${OUTPUTDIR}/${PRESET}
cp -r $PREFIX.bam ${OUTPUTDIR}/${PRESET}
cp -r $PREFIX.bam.bai ${OUTPUTDIR}/${PRESET}
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
