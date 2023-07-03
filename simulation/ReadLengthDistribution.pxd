cdef class RLDfactory_gamma:
    cdef float __alpha
    cdef float __loc
    cdef float __beta
    cdef int nextl(self)


cdef RLDfactory_gamma get_rld_factory(int mean, float alpha, float loc, float beta, str rldfile)