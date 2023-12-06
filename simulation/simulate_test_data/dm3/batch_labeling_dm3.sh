#! /bin/bash

[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/benchmark/label_filelist
REPEAT='/data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.rmsk_200.bed'
GAP='/data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3_gap.bed'
TE='/data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa'
BLACKLIST='/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/For_ML/germ_v4/Dedup/BlackList.bed'
TEMP3='/data/tusers/zhongrenhu/for_SMS/bin/demo3/develop'
germModel='/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/For_ML/germ_v4/Dedup/dm3_G_1V1/'
somaModel='/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/For_ML/soma_v4/Dedup/dm3_S_1V30/'


germlineList=(1 3 5)
for i in ${germlineList[*]}
do
  str=`sed "${i}q;d" $FILELIST`
  toks=($str)
  BAM=${toks[0]}
  WORKDIR=`dirname $BAM`
  OUTPUTDIR=${WORKDIR}/result_no_secondary
  WORKDIR=`dirname $WORKDIR`
  [ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR
  cd $OUTPUTDIR && echo "Changing wd: " $(pwd) && echo "Labeling for ${BAM}"

  /data/tusers/zhongrenhu/for_SMS/bin/simulation/label/label_protocol.sh \
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
    --subSize 5 \
    --germModel $germModel \
    --somaModel $somaModel

  echo "Done." && echo ""
done


somaticlineList=(2 4 6)
for i in ${somaticlineList[*]}
do
  str=`sed "${i}q;d" $FILELIST`
  toks=($str)
  BAM=${toks[0]}
  WORKDIR=`dirname $BAM`
  OUTPUTDIR=${WORKDIR}/result_no_secondary
  WORKDIR=`dirname $WORKDIR`
  [ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR
  cd $OUTPUTDIR && echo "Changing wd: " $(pwd) && echo "Labeling for ${BAM}"

  /data/tusers/zhongrenhu/for_SMS/bin/simulation/label/label_protocol.sh \
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

  echo "Done." && echo ""
done
