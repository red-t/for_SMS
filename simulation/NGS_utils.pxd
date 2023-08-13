from fastaIO cimport readAllTuples, rc, get_length_list
from fastqIO cimport FastqPairWriter
from Mutator cimport NGS_Mutator
from CoverageGenerator cimport RandomReads

cpdef generate_PE(str pop_gen,
                  float error_rate,
                  int nread,
                  int insertsize,
                  int stddev,
                  int rlen,
                  str fq1,
                  str fq2,
                  float chimera)