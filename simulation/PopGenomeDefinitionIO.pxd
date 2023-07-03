from TESequenceBuilder cimport SequenceContainer

cdef class PopGenDefinitionReader:
    cdef public:
        str __filename
        str __chasis
        int insertions
        int popsize
        list tuples
        SequenceContainer __sc

    cdef __parse_chasis(self, str line)
    cpdef str  get_chasis(self)
    cdef list read_tuples(self)
    cpdef list read_transposed(self)