from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = [
    Extension(name = "CoverageGenerator", sources = ["CoverageGenerator.pyx"]),
    Extension(name = "define_pgdf_utils", sources = ["define_pgdf_utils.pyx"]),
    Extension(name = "fastaIO", sources = ["fastaIO.pyx"]),
    Extension(name = "fastqIO", sources = ["fastqIO.pyx"]),
    Extension(name = "Mutator", sources = ["Mutator.pyx"]),
    Extension(name = "NGS_utils", sources = ["NGS_utils.pyx"]),
    Extension(name = "PopGenomeDefinitionIO", sources = ["PopGenomeDefinitionIO.pyx"]),
    Extension(name = "ReadLengthDistribution", sources = ["ReadLengthDistribution.pyx"]),
    Extension(name = "TEInsert", sources = ["TEInsert.pyx"]),
    Extension(name = "TESequenceBuilder", sources = ["TESequenceBuilder.pyx"]),
    Extension(name = "TGS_utils", sources = ["TGS_utils.pyx"]),
    Extension(name = "BPG_utils", sources = ["BPG_utils.pyx"])
]

setup(ext_modules=cythonize(ext, language_level=3))