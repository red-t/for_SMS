#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH -c 24
#SBATCH --array=2,4,6
#SBATCH --partition=12hours
#SBATCH --output=./logs/Label/GRCh38/label-log-%A-%a.out


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
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/benchmark/label_filelist
[ -z $CONTIGS ] && CONTIGS=(`cut -f 1 /data/tusers//zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_germ/HG02716_template.fa.fai`)

REPEAT='/data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/rmsk_200.bed'
GAP='/data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/gap.bed'
TE='/data/tusers/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38.transposon.fa'
BLACKLIST='/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/For_ML/germ_v4/Dedup/BlackList.bed'
TEMP3='/data/tusers/zhongrenhu/for_SMS/bin/demo3/develop'
germModel='/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/For_ML/germ_v4/Dedup/GRCh38_G_1V1/'
somaModel='/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/For_ML/soma_v4/Dedup/GRCh38_S_1V30/'

### Get parameters from FILELIST ###
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
BAM=${toks[0]}
WORKDIR=`dirname $BAM`
OUTPUTDIR=${WORKDIR}/result_no_secondary
WORKDIR=`dirname $WORKDIR`


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
source /zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh
conda activate TEMP3


########################
### Process the data ###
########################
echo "Labeling for ${BAM}"
bash /data/tusers/zhongrenhu/for_SMS/bin/simulation/label/label_protocol.sh \
    -b $BAM \
    -d $WORKDIR \
    -r $REPEAT \
    -g $GAP \
    -P 10 \
    -T 5 \
    --refTe $TE \
    --blacklist $BLACKLIST \
    --minLen 100 \
    --TEMP3 $TEMP3 \
    --subSize 50 \
    --germModel $germModel \
    --somaModel $somaModel
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
