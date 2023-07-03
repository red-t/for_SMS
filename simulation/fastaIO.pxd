cpdef str rc(str seq)

cpdef tuple load_chasis(str fastafile)

cpdef list get_length_list(str inputFile)

cpdef list readAllTuples(str fastafile)

cdef class FastaWriter:
    cdef str __filename
    cdef object __filehandle
    cdef int __seqleng

    cpdef write(self, str name, str seq)