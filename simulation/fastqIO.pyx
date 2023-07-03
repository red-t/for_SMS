cdef class FastqWriter:
    def __init__(self, str file):
        self.__filename = file
        self.__filehandle = open(file,"w")

    cdef write(self, str header, str seq):
        fh = self.__filehandle
        fh.write("@" + header + "\n")
        fh.write(seq + "\n")
        fh.write("+" + header + "\n")
        fh.write("I" * len(seq) + "\n")
    
    def close(self):
        self.__filehandle.close()


cdef class FastqPairWriter:
    def __init__(self, str file1, str file2):
        self.__fqw1 = FastqWriter(file1)
        self.__fqw2 = FastqWriter(file2)
    
    cdef write(self, str header, str seq1, str seq2):
        self.__fqw1.write(header, seq1)
        self.__fqw2.write(header, seq2)

    def close(self):
        self.__fqw1.close()
        self.__fqw2.close()