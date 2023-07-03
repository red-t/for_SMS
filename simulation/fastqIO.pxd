cdef class FastqWriter:
    cdef str __filename
    cdef object __filehandle
    cdef write(self, str header, str seq)


cdef class FastqPairWriter:
    cdef FastqWriter __fqw1
    cdef FastqWriter __fqw2
    cdef write(self, str header, str seq1, str seq2)