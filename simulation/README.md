## 1. Environment Configuration

### 1.1 Install dependencies

```shell
# Create enviroment
mamba create -n simulation
mamba activate simulation

# VISOR
mamba install VISOR

# bcftools
mamba install bcftools

# samtools
mamba install samtools

# scipy
mamba install scipy

# pysam
mamba install pysam

# cython
pip install cython
```

### 1.2 Compile cython scripts

```shell
# Your current working directory should be <path/to/simulation>
python setup.py build_ext -i && rm -r build && rm *c
```

## 2. Build template genome

Based on **reference genome** and **polymorphisms resources** (from large population scale research, like [DGRP2](http://dgrp2.gnets.ncsu.edu/data.html), [GGVP](https://www.internationalgenome.org/data-portal/data-collection/ggvp-grch38), [1000 genomes project phase 3](https://www.internationalgenome.org/category/phase-3/)), this step aims to build template haplotype genome with SNV correspond to specific strain/sample.

### 2.1 For Human

```shell
# Using GRCh38 as reference genome
# Using SNVs of HG02716 from GGVP
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_germ

./build_template_genome.sh -n HG02716 -s F -f /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_X.fa -m /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GRCh38_no_alt_Y.fa -v /data/tusers.ds/zhongrenhu/for_SMS/reference/GRCh38.p13/GGVP/
```

The main results under `/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/template/HG02716_germ`:

```shell
# Only use HG02716_template_h1.fa later

HG02716_template_h1.fa			# simulated haplotype 1 template genome
HG02716_template_h1.fa.fai
HG02716_template_h2.fa			# simulated haplotype 2 template genome
HG02716_template_h2.fa.fai
```

### 2.2 For Drosophila

Run script `build_template_genome_unphased.sh` to build template genome, for example, with dm3:

```shell
# Using dm3 as reference genome
# Using SNVs of line_28 from DGRP
# Current working directory is /data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/template/line_28

./build_template_genome_unphased.sh -n line_28 -r /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.fa -v /data/tusers.ds/zhongrenhu/for_SMS/reference/dm3/DGRP/
```

The main results under `/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm3/template/line_28`:

```shell
line_28_template.fa        # simulated template genome
line_28_template.fa.fai
```

## 3. Build population genome and NGS, TGS data

Based on **"template genome"**, this step aims to build **"population genome"** for each chromosome and generate raw NGS and TGS reads from it.

### 3.1 Running on local machine

You can finish this step on local machine by running `simulation_protocol_local.sh`:

```shell
./simulation_protocol_local.sh -d ./ -r line_28_template.fa -t /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.transposon_for_simulaTE.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 50 --min-distance 9000 --depth 50 --species fly --protocol ccs
```

The main results under working directory:

```shell
NGS_*.fastq					# simulated NGS raw data
TGS.fasta					# simulated TGS raw data
chr2L/chr2L.ins.summary     # summary information of each simulated insertion
chr2L/chr2L.ins.sequence    # sequence of each simulated insertion
```

### 3.2 Running on blades

Actually, I fnished this step via blades, see `exmaples`.

## 4. Details of output

Each row in `chr2L.ins.summary` represents a simulated insertion:

```shell
# Column	ColumnID	Description

# 1			chrom		chromosome
# 2			start		left breakpoint of the insertion (1-based)
# 3			end			right breakpoint of the insertion (1-based)
# 4			insId		insertion id
# 5			insPgd		insertion "population genome definition", which representing the structure of inserted sequence, see https://sourceforge.net/p/simulates/wiki/describing_TE_sequences
# 6			strand		orientation of the inserted sequence (related to TE CSS)
# 7			frequency	insertion frequency in the "population genome"
# 8			tsdStart	start of TSD (1-based)
# 9			tsdEnd		end of TSD (1-based)
```

Each row in `chr2L.ins.sequence` represents the sequence of a simulated insertion:

```shell
# Column	ColumnID	Description

# 1			insId		insertion id, same as chr2L.ins.summary
# 2			insSeq		insertion sequence, including flanking region
# 							1. the structure is 'flank -- TSD -- insert -- TSD -- flank'
# 							2. if (insertion_start < 2000), l=start, or l=2000. And TSD is insSeq[l:l+len(tsd)]
# 							3. if (insertion_end + 2000 >= len(ref)+1), r=(len(ref) - 1 - insertion_end), or r=2000. And TSD is insSeq[-r-len(tsd):-r]
```
