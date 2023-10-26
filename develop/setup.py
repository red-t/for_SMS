from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext = [
    Extension(name = "TEMP3.htslib_external",
              sources = ["TEMP3/htslib_external.pyx"],
              libraries = ["hts"]),

    Extension(name = "TEMP3.AlignmentFileIO",
              sources = ["TEMP3/AlignmentFileIO.pyx"],
              libraries = ["hts"]),

    Extension(name = "TEMP3.Cluster",
              sources = ["TEMP3/Cluster.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c"],
              libraries = ["hts"]),

    Extension(name = "TEMP3.ParallelTemplate",
              sources = ["TEMP3/ParallelTemplate.pyx",  "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c"],
              libraries = ["hts"])
    ]

setup(ext_modules=cythonize(ext, language_level=3))