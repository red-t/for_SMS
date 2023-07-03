from fastaIO cimport readAllTuples, rc, get_length_list, FastaWriter
from Mutator cimport PacBioMutator
from CoverageGenerator cimport RandomReads
from ReadLengthDistribution cimport RLDfactory_gamma, get_rld_factory

cpdef generate_TGS(str pop_gen,
                   float error_rate,
                   str err_frac,
                   int nread,
                   int rlen,
                   float alpha,
                   float loc,
                   float beta,
                   str rldfile,
                   str outfa,
                   int tgs_minl,
                   int tgs_maxl)