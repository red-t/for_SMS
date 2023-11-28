----

## Step 1: build unphased template genome

Based on **reference genome** and **polymorphisms resources** (from large population scale research, like [DGRP2](http://dgrp2.gnets.ncsu.edu/data.html), [GGVP](https://www.internationalgenome.org/data-portal/data-collection/ggvp-grch38), [1000 genomes project phase 3](https://www.internationalgenome.org/category/phase-3/)), this step aims to build template haplotype genome with SNV correspond to specific strain/sample.

Make sure you have the paths of the executable program `VISOR` and `bcftools` in your environment variable `PATH`:
```
# VISOR and bcftools can be installed by conda
conda create -n VISOR
conda install VISOR -c conda-forge
conda install bcftools -c conda-forge
```

#### For Human:
Run script `build_template_genome.sh` to build template genome, for example, with GRCh38:
```
# HG02716 是 GGVP 当中的女性个体之一
build_template_genome.sh -n HG02716 -s F -f /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_X.fa -m /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_Y.fa -v /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GGVP/
```

And the main results look like:
```
# 路径 /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/HG02716_templateswithsnp
HG02716_template_h1.fa			# simulated haplotype 1 template genome
HG02716_template_h1.fa.fai		# index of it
HG02716_template_h2.fa			# simulated haplotype 2 template genome
HG02716_template_h2.fa.fai		# index of it
```
(后续可能只需要使用 female 的 haplotype 1 作为 template genome)

#### For Fruitfly:
Run script `build_template_genome_unphased.sh` to build template genome, for example, with dm3:
```
# line_28 是 DGRP 当中的品系之一
build_template_genome_unphased.sh -n line_28 -r /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.fa -v /data/tusers.ds/zhongrenhu/for_SMS/reference/dm3/DGRP/
```

And the main results look like:
```
# 路径 /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28
line_28_template.fa        # simulated haplotype template genome
line_28_template.fa.fai    # the index of it
```

----

## Step 2: simulate population genome and NGS, TGS data

Based on the **"template genome"** build at step1, here we will build **"population genome"** for each chromosomes and generate NGS and TGS from it.

Make sure you have these dependencies:
```
1. samtools
conda install samtools -c conda-forge

2. scipy
conda install scipy -c conda-forge

3. pysam
conda install pysma -c conda-forge

4. cython
pip install cython
```

Then compiling
```
# make sure you are now in path/to/simulation
python setup.py build_ext -i && rm -r build && rm *c
```

Then run script `simulation_protocol.sh`:

#### For Fruitfly
```
# 此时的工作路径是 /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/line_28
simulation_protocol.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 10000 -R 0 --sub-N 50 --germline-count 500 --avg-somatic-count 20 --min-distance 500 --depth 50
```

And the main results look like:
```
line_28_*.fastq             # simulated NGS raw data
line_28_pacbio.fasta        # simulated TGS raw data
chr2L/chr2L.ins.summary     # summary information of each simulated insertion
chr2L/chr2L.ins.sequence    # sequence of each simulated insertion
```

----

## More details of outputs

***chr2L.ins.summary*** 记录每个 simulated  insertion 的信息，每一行的内容如下:
```
# chrom l_bp r_bp ins_id ins_pgd orient frequency tsd_start tsd_end
chr2L   19728056        19728057        33009~FBgn0005384_3S18~UF       $15+6bp +       0.0001  19728051        19728056

- l_bp
	Insertion 在 template genome 当中的 "left breakpoint", 1-based index.

- r_bp
	Insertion 在 template genome 当中的 "right breakpoint", 1-based index. Insertion 插入的位置就是 l_bp 和 r_bp 之间.

- ins_id
	Insertion 的 id, 其中的 'U / T' 表示 'un- / truncated', 而其中的 'F / R' 表示 'forward / reversed'.

- ins_pgd
	Insertion 的 "population genome definition", 保留了插入片段的结构, 解读方式参考: https://sourceforge.net/p/simulates/wiki/describing_TE_sequences

- orient
	插入片段的方向 (相对于 transposon CSS 而言, 只表示最外层 insertion 的方向).

- frequency
	Insertion 的 frequency, 是指该 insertion 在 population genome 当中的 frequency.

- tsd_start
	TSD 在 template genome 当中的起始位点, 1-based index.

- tsd_end
	TSD 在 template genome 当中的终止位点, 1-based index. 通过 samtools faidx template.fa chrXX:tsd_start-tsd_end 可以获取 TSD 序列
```

***chr2L.ins.sequence*** 记录每个 simulated  insertion 的序列信息，每一行的内容如下:
```
# ins_id seq
chr2L_hg458;33007~FBgn0003055_P_element~UR ATGTGAATACGAAATT...

- ins_id
	Insertion 的 id, ';' 符号之前的内容表示这条序列的来源

- seq
	Insertion 的序列, 包含两端的 flank region, 结构为 'flank -- TSD -- insert -- TSD -- flank'
	如果该 insertion 的 l_bp < 2000, l=l_bp, 否则 l=2000, 相应的: seq[l:l+len(tsd)] 为 TSD 序列
	如果该 insertion 的 r_bp + 2000 >= len(ref)+1, r=len(ref)-1-l_bp, 否则 r=2000, 相应的: seq[-r-len(tsd):-r] 为 TSD 序列
```

