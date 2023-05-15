from distutils.core import setup, Extension
from Cython.Build import cythonize
import pysam

## V1 ##
# ext = [
#     Extension(name="htslib_external", sources=["htslib_external.pyx"],
#                 include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
#     Extension(name="SegmentParser", sources=["SegmentParser.pyx"],
#                 include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
#     Extension(name="AlignmentFileIO", sources=["AlignmentFileIO.pyx"],
#                 include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
#     Extension(name="Cluster", sources=["Cluster.pyx"],
#                 include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
#     Extension(name="ParallelTemplate", sources=["ParallelTemplate.pyx"],
#                 include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"])
#     ]

## V2 ##
# ext = [
#     Extension(name="htslib_external", sources=["htslib_external.pyx"],
#                 include_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/include"]+pysam.get_include(), libraries=["hts"], library_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/lib"]),
#     Extension(name="SegmentParser", sources=["SegmentParser.pyx"],
#                 include_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/include"]+pysam.get_include(), libraries=["hts"], library_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/lib"]),
#     Extension(name="AlignmentFileIO", sources=["AlignmentFileIO.pyx"],
#                 include_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/include"]+pysam.get_include(), libraries=["hts"], library_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/lib"]),
#     Extension(name="Cluster", sources=["Cluster.pyx"],
#                 include_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/include"]+pysam.get_include(), libraries=["hts"], library_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/lib"]),
#     Extension(name="ParallelTemplate", sources=["ParallelTemplate.pyx"],
#                 include_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/include"]+pysam.get_include(), libraries=["hts"], library_dirs=["/Users/hzr/opt/anaconda3/envs/pysam/lib"])
#     ]

## V3 ##
# incl_dirs 当中包含 '/Users/hzr/opt/anaconda3/envs/pysam/lib/python3.9/site-packages/pysam/include/htslib'
# 该路径与 /Users/hzr/opt/anaconda3/envs/pysam/include 一样，都包含更深一层的 htslib 子目录(其中包含 htslib 的一些头文件)
# 前者是通过 conda 安装 htslib 生成的；后者是通过 conda 安装 pysam 生成的，它们的版本一样
# 因此，这里直接通过 pysam.get_include() 获取 include_dirs 也是可以的，没必要重复。(当然，如果版本不同的话就不行了)
incl_dirs = pysam.get_include()
lib_dirs = ["/Users/hzr/opt/anaconda3/envs/pysam/lib"]

ext = [
    Extension(name = "htslib_external",
              sources = ["htslib_external.pyx"],
              include_dirs = incl_dirs,
              libraries = ["hts"],
              library_dirs = lib_dirs),

    Extension(name = "SegmentParser",
              sources = ["SegmentParser.pyx"],
              include_dirs = incl_dirs,
              libraries = ["hts"],
              library_dirs = lib_dirs),

    Extension(name = "AlignmentFileIO",
              sources = ["AlignmentFileIO.pyx"],
              include_dirs = incl_dirs,
              libraries = ["hts"],
              library_dirs = lib_dirs),

    Extension(name = "Cluster",
              sources = ["Cluster.pyx"],
              include_dirs = incl_dirs,
              libraries = ["hts"],
              library_dirs = lib_dirs),

    Extension(name = "ParallelTemplate",
              sources = ["ParallelTemplate.pyx"],
              include_dirs = incl_dirs,
              libraries = ["hts"],
              library_dirs = lib_dirs)
    ]

setup(ext_modules=cythonize(ext, language_level=3))  # 指定Python3