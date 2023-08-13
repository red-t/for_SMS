from fastaIO cimport readAllTuples, rc, get_length_list, FastaWriter
from Mutator cimport TGS_Mutator
from CoverageGenerator cimport RandomReads
from ReadLengthDistribution cimport RLDfactory

cpdef generate_PACBIO(str pop_gen,
                      int nread,
                      int minl,
                      int maxl,
                      str outfa,
                      str protocol)


cpdef generate_ONT(str pop_gen,
                   int nread,
                   int minl,
                   int maxl,
                   str outfa,
                   str protocol)