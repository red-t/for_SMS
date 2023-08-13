cdef class TGS_Mutator:
    cdef float   __er
    cdef dict    __tr
    cdef list    __ins
    cdef float   __misf
    cdef float   __delf

    cdef str mutateseq(self, str seq)


cdef class NGS_Mutator:
    cdef float   __er
    cdef dict    __tr
    
    cdef int __getPoisson(self, float lam)
    cdef str mutateseq(self, str seq)


cdef class ExhaustiveSeqMutator:
    cdef float   __er
    cdef dict    __tr
    
    cdef str mutateseq(self, str seq)