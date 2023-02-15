----
#### Dependencies

**1. pysam 0.20.0**
```
# install with conda (recommended)
conda config --add channels r
conda config --add channels bioconda
conda install pysam=0.20.0

# install through pypi
pip install
```

**2. htslib 1.16**
```
# download
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -jxvf htslib-1.16.tar.bz2 && mv htslib-1.16 htslib

# compile (at present, don't need to install)
cd htslib
./configure --prefix=/xx/xx/xx
make

# or install with conda (not support yet)
conda install htslib=1.16
```

**3. cython**
```
# install through pypi
pip install Cython --install-option="--no-cython-compile"

# maybe conda is also ok?
```

----
#### Compile

```
# now your workdir should be path/to/TEMP3
python setup.py build
mv build/lib*/* . && rm -r build/ && rm *c

# once compiled, only the executable file *.so is needed.
```

----
#### Usage

```
# now your workdir should be path/to/TEMP3
import AlignmentFileIO

fn = "ex2.bam" # should be indexed
nthread = 5
stid = 0
maxtid = 15
bf = AlignmentFileIO.BamFile(fn, nthread)
a = []
for b in bf.fetch(stid, maxtid):
	if len(b):
		a.extend(b)

len(a)
bf.close()
```