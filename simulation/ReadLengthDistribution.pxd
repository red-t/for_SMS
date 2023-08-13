cdef class RLDfactory:
    cdef float __alpha
    cdef float __loc
    cdef float __scale
    cdef int gamma_nextl(self)
    cdef int expon_nextl(self)