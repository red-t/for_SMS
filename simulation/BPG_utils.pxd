from fastaIO cimport readAllTuples, load_chasis, FastaWriter
from PopGenomeDefinitionIO cimport PopGenDefinitionReader
from TESequenceBuilder cimport SequenceContainer
from TEInsert cimport insertSequences


cpdef build_popg(str te_fasta,
                 str pgdf,
                 str ref_fasta,
                 int sub_idx,
                 int sub_size,
                 str output,
                 str outseqf)