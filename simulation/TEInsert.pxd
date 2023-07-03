cdef class TESequence:
    cdef str sequence
    cdef int id
    cdef int tsd


cpdef str insertSequences(str ref, list posinstuples, str fout_name=*, list seqidtups=*)