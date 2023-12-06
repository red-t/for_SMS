#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH -c 24
#SBATCH --array=1-50%10
#SBATCH --partition=4hours
#SBATCH --output=./logs/label-log-%A-%a.out

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
[ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/dm3/line_21_0/label_filelist
# [ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/HG02716_S0/label_filelist
# [ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/HG02716_G0/label_filelist

[ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/dm3/template/line_21/line_21_template.fa.fai`)
# [ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_soma/HG02716_template.fa.fai`)
# [ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers.ds/zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_germ/HG02716_template.fa.fai`)

### Get parameters from FILELIST ###
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
BAM=${toks[0]}
RMSK=${toks[1]}
GAP=${toks[2]}
TE=${toks[3]}
WORKDIR=`dirname $BAM`
OUTPUTDIR=${WORKDIR}/result_no_secondary


###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR

### Copy inputs to temp dir ###
echo "Copying input file ${BAM} to temp directory..."
for contig in ${CONTIGS[*]}; do
    cp -r $WORKDIR/$contig .
done
cp -r $OUTPUTDIR .
cd result_no_secondary && echo "Changing wd: " $(pwd)
echo "Done."
echo ""


##################################
### Activate conda environment ###
##################################
source /data/tusers/zhongrenhu/Software/anaconda3/etc/profile.d/conda.sh
conda activate TEMP3


########################
### Process the data ###
########################
echo "Labeling for ${BAM}"
bash /data/tusers.ds/zhongrenhu/for_SMS/bin/demo3/test_with_simulation/label_protocol.sh -b $BAM -d $TMPDIR -r $RMSK -g $GAP -T $TE -p 10 -t 5 -l 100 -s 50
# bash /data/tusers.ds/zhongrenhu/for_SMS/bin/demo3/test_with_simulation/label_protocol.sh -b $BAM -d $TMPDIR -r $RMSK -g $GAP -T $TE -p 10 -t 5 -l 100 -s 20 -i 19
# bash /data/tusers.ds/zhongrenhu/for_SMS/bin/demo3/test_with_simulation/label_protocol.sh -b $BAM -d $TMPDIR -r $RMSK -g $GAP -T $TE -p 10 -t 5 -l 100 -s 5
echo ""


################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp AllIns* $OUTPUTDIR
cp *bed $OUTPUTDIR
cp *txt $OUTPUTDIR
cp merge* $OUTPUTDIR
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
