from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

include_dirs = ["./AIList/", numpy.get_include()]
ext = [
    Extension('ailist',
              sources = ['./AIList/ail.pyx', './AIList/AIList.c'],
              depends = ['AIList.h', 'khash.h', 'kseq.h', 'ail.pyx'], 
              include_dirs = include_dirs)
]

setup(ext_modules=cythonize(ext, language_level=3))
