from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext = [
    Extension(name = "TEMP3.FileIO",
              sources = ["TEMP3/FileIO.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c", "TEMP3/src/io_utils.c"],
              libraries = ["hts"]),

    Extension(name = "TEMP3.Cluster",
              sources = ["TEMP3/Cluster.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c", "TEMP3/src/io_utils.c"],
              libraries = ["hts"]),
    
    Extension(name = "TEMP3.Assemble",
              sources = ["TEMP3/Assemble.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c", "TEMP3/src/io_utils.c"],
              libraries = ["hts"]),
    
    Extension(name = "TEMP3.Annotate",
              sources = ["TEMP3/Annotate.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c", "TEMP3/src/io_utils.c", "TEMP3/src/anno_utils.c", "TEMP3/src/post_filter.c"],
              libraries = ["hts"]),
    
    Extension(name = "TEMP3.ParallelModule",
              sources = ["TEMP3/ParallelModule.pyx", "TEMP3/src/AIList.c", "TEMP3/src/seg_utils.c", "TEMP3/src/cluster_utils.c", "TEMP3/src/io_utils.c", "TEMP3/src/anno_utils.c", "TEMP3/src/post_filter.c", "TEMP3/src/ltr_utils.c"],
              libraries = ["hts"]),
    ]

setup(ext_modules=cythonize(ext, language_level=3))