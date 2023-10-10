----
#### Dependency
**1. Autogluon**
```
conda create -n TEMP3 python=3.10

conda activate TEMP3

conda install -c conda-forge mamba

mamba install -c conda-forge autogluon
```

**2. Samtools (1.7)**
```
# 最新版本的应该就是 1.7

conda install samtools -c conda-forge
```

**3. Minimap2 (2.26-r1175)**
```
# 最新版本的应该就是 2.26-r1175

conda install minimap2 -c conda-forge
```

**4. Cython (0.29.33)**
```
# 最新版本的应该就是 0.29.33

pip install cython
```

----
#### Compile
```
# Now your workdir should be path/to/develop
# Once compiled, only the executable file *.so is needed.

python setup.py build_ext -i
rm -r build/ && rm TEMP3/*c
```

----
#### Usage
```
python path/to/develop/TEMP3.py -b BAM -r REPEAT.bed -g GAP.bed -T TE.fa -p NUM_PROCESS -t NUM_THREAD -l 100
```

---
#### Result
```
# Now, each column of the *clt.txt is:
# Column 1 --- chrom
# Column 2 --- start
# Column 3 --- end
# Column 4 --- cluster id
# Column 5 --- number of segment (normalized by background depth)
# Column 6 --- strand
# Column 7 --- start index (index of corresponding segments)
# Column 8 --- end index
# Column 9 --- nseg (normalized by background depth)
# Column 10 --- strand flag (the orientation of this insertion)
					1:forward
					2:reverse
					other:unkown
# Column 11 --- single flag (whether the cluster only have single support read)
					0:multiple support reads
					1: single read with 1 alignment
					2: single read with 2 alignments
# Column 12 --- cluster location flag
					1: at normal region
					2: at repeat/gap boundary
					4: inside repeat/gap
# Column 13 --- number of segment type
# Column 14 --- entropy
# Column 15 --- balance ratio
# Column 16 --- fraction of segments with low mapq (<5)
# Column 17 --- fraction of "dual-clip" alignments
# Column 18 --- fraction of segments with loc_flag=1
# Column 19 --- fraction of segments with loc_flag=2
# Column 20 --- fraction of segments with loc_flag=4
# Column 21 --- fraction of segments with loc_flag=8
# Column 22 --- fraction of segments with loc_flag=16
# Column 23 --- average mapq of this cluster
# Column 24 --- average per base alignment score
# Column 25 --- average query aligned fraction
# Column 26 --- average per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
# Column 27 --- average per-base gap-compressed divergence (normalized background)
# Column 28 --- background div
# Column 29 --- background de
# Column 30 --- background depth
# Column 31 --- background read length
# Column 32 --- fraction of TE-aligned segment
# Column 33 --- flag (if the cluster pass TE-alnfrac threshold)
					0: passed
					1: low alnfrac
# Column 34 --- TE (tid of corresponding TE)
```