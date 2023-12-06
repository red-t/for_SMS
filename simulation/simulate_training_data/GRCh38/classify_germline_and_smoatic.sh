#! /bin/bash

[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/v4/HG02716_0_ccs/label_filelist
[ -z $MIN_FREQ ] && MIN_FREQ=0.01

for i in $(seq 1 50) $(seq 251 300) $(seq 401 450)
do
  str=`sed "${i}q;d" $FILELIST`
  toks=($str)
  BAM=${toks[0]}
  WORKDIR=`dirname $BAM`
  OUTPUTDIR=${WORKDIR}/result_no_secondary
  WORKDIR=`dirname $WORKDIR`
  cd $OUTPUTDIR
  cat $WORKDIR/*/*summary > tmp.summary

  python /data/tusers/zhongrenhu/for_SMS/bin/demo3/v4_test/classify_germline_and_smoatic.py \
    --summary tmp.summary \
    --input TP.bed \
    --output tmp_TP.bed \
    --freq ${MIN_FREQ}
  
  cut -f 1-4 TP_clt.txt > t1
  cut -f 5 tmp_TP.bed > t2
  cut -f 6-33 TP_clt.txt > t3
  paste t1 t2 t3 > tmp_TP_clt.txt
  awk -v freq=${MIN_FREQ} '$5>freq{print $0}' tmp_TP_clt.txt > TP_clt_G.txt
  awk -v freq=${MIN_FREQ} '$5==freq{print $0}' tmp_TP_clt.txt > TP_clt_S.txt
  rm tmp.summary && rm tmp_TP.bed && rm tmp_TP_clt.txt && rm t1 && rm t2 && rm t3
done
