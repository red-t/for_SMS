----
#### Dependency
**1. Autogluon**
```
conda create -n TEMP3 python=3.10.13

conda activate TEMP3

conda install -c conda-forge mamba

mamba install -c conda-forge autogluon
```

**2. Samtools (1.17)**
```
conda install samtools=1.17 -c conda-forge

# or
mamba install samtools=1.17 -c conda-forge
```

**3. Minimap2 (2.26-r1175)**
```
# 最新版本的应该就是 2.26-r1175

conda install minimap2 -c conda-forge

# or
mamba install minimap2 -c conda-forge
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
python TEMP3.py -b BAM -r REPEAT.bed -g GAP.bed -B BlackList.bed -T TE.fa --germ GERM_MODEL --soma SOMA_MODEL -p NUM_PROCESS -t NUM_THREAD
```

---
#### Result
```
# Now, each column of the *clt.txt is:

Column	Value	Description

1	chrom				chromosome
2	refStart			cluster start on reference sequence (0-based, included)
3	refEnd				cluster end on reference sequence (0-based, not-included)
4	cltID				cluster ID
5	numSeg				number of segments in the cluster (normalized by bg depth)
6	strand				cluster orientation
7	startIndex			start index in segments array (0-based, include)
8	endIndex			end index in segments array (0-based, not-include)
9	numSeg				number of segments in the cluster (normalized by bg depth)
10	directionFlag		bitwise flag representing cluster direction
							1: forward
							2: reverse
							other: unkown
11	cltType				cluster type
							0: germline (multiple support reads)
							1: somatic (1 support read & 1 alignment)
							2: somatic (1 support read & 2 alignments)
12	locationType		bitwise flag representing cluster location
							1: inside normal region
							2: at repeat/gap boundary
							4: inside repeat/gap
13	numSegType			number of different segment types
14	entropy				entropy based on fraction of different type segments
15	balanceRatio		balance ratio based on number of left- & right-clip segments
16	lowMapQualFrac		fraction of segments with low mapQual (<5)
17	dualClipFrac		fraction of "dual-clip" alignments
18	alnFrac1			fraction of segments with alnLocationType=1
19	alnFrac2			fraction of segments with alnLocationType=2
20	alnFrac4			fraction of segments with alnLocationType=4
21	alnFrac8			fraction of segments with alnLocationType=8
22	alnFrac16			fraction of segments with alnLocationType=16
23	meanMapQual			mean mapQual of cluster
24	meanAlnScore		mean per-base alignment score (based on teAlignments)
25	meanQueryMapFrac	mean query mapped fraction (based on teAlignments)
26	meanDivergence		mean per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
27	bgDiv				background divergence (for normalization)
28	bgDepth				background depth (for normalization)
29	bgReadLen			background read length
30	teAlignedFrac		fraction of TE-aligned segments
31	teTid				majority TE-tid of cluster
32	isInBlacklist		whether cluster intersects with blacklist
33	probability			the probability of the cluster to be a positive insertion
```