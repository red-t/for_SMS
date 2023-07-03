cdef class RandomReads:
    cdef list __reads
    cdef list getreadtable(self, int readnum, list pgll)
    cdef int get_reads(self, int index)