########################################################
## global variables
# maximum genomic coordinace
# for some reason, using 'int' causes overflow
cdef int MAX_POS = (1 << 31) - 1
cdef dict SEG_DICT = {}