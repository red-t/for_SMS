import pysam

cpdef str rc(str seq):
    cdef:
        dict complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N',
                           'a': 'T', 'c': 'G', 'g': 'C', 't': 'A', 'n':'N',}
        str base, ret
    
    ret = ''.join([complement[base] for base in seq[::-1]])
    return ret


cpdef tuple load_chasis(str fastafile):
    cdef str name, seq
    with pysam.FastxFile(fastafile) as fh:
        for entry in fh:
            name = entry.name
            seq = entry.sequence
    
    return (name, seq)


cpdef list get_length_list(str inputFile):
    cdef:
        list ll = []

    with pysam.FastxFile(inputFile) as fh:
        for entry in fh:
            ll.append(len(entry.sequence))
    
    return ll


cpdef list readAllTuples(str fastafile):
    cdef:
        list ft = []
        str name, seq
    
    with pysam.FastxFile(fastafile) as fh:
        for entry in fh:
            ft.append((entry.name, entry.sequence))
    
    return ft


cdef class FastaWriter:
    """
    Write the content to a fasta file
    """
    def __init__(self, str file, int seqleng):
        self.__filename = file
        self.__filehandle = open(file,"w")
        self.__seqleng = seqleng
    
    cpdef write(self, str name, str seq):
        cdef:
            int sl = self.__seqleng
            int c = 0
        self.__filehandle.write(">" + name + "\n")
        while(c < len(seq)):
            self.__filehandle.write(seq[c:c+sl] + "\n")
            c += sl
    
    def close(self):
        self.__filehandle.close()