from distutils.core import setup, Extension
from Cython.Build import cythonize
import pysam

ext = [
    Extension(name="htslib_external", sources=["htslib_external.pyx"],
                include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
    Extension(name="SegmentParser", sources=["SegmentParser.pyx"],
                include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"]),
    Extension(name="AlignmentFileIO", sources=["AlignmentFileIO.pyx"],
                include_dirs=["./htslib"]+pysam.get_include(), libraries=["hts"], library_dirs=["./htslib"])
    ]

setup(ext_modules=cythonize(ext, language_level=3))  # 指定Python3