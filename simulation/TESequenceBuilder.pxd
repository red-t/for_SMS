from fastaIO cimport rc
from TEInsert cimport TESequence, insertSequences
from Mutator cimport ExhaustiveSeqMutator

cdef class SequenceContainer:
    cdef list __ssl
    cdef dict __sc

    cdef addDefinition(self, str definition)
    cpdef TESequence getTESequence(self, str id)
    cdef TESequence __process_nested(self, TESequence teseq, str nestdef)
    cdef list __recursionsafe_split(self, str nestdef)
    cdef TESequence __process_right(self, TESequence teseq, str right)
    cdef TESequence __fixstrand(self, TESequence teseq, str strand)
    cdef TESequence __process_deletions(self, TESequence teseq, str deldef)
    cdef TESequence __commence_deleting(self, TESequence teseq, set delset)
    cdef TESequence __teseqfromHash(self, str seqid)
    cdef int __getPosition(self, str pos, TESequence teseq)