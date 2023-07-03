import random

cdef class RandomReads:
    def __init__(self, int readnum, list pgll):
        self.__reads = self.getreadtable(readnum, pgll)

    cdef list getreadtable(self, int readnum, list pgll):
        cdef:
            list table = [0,] * len(pgll)
            int i, r, n = len(pgll)-1
        for i in range(0, readnum):
            r = random.randint(0, n)
            table[r] += 1

        return table

    cdef int get_reads(self, int index):
        """
        zero based index of the genome
        """
        return self.__reads[index]