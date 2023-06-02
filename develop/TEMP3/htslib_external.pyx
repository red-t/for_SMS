########################
### global variables ###
########################

cdef int MAX_POS = (1 << 31) - 1
cdef str TEXT_ENCODING = 'utf-8'
cdef str ERROR_HANDLER = 'strict'



# cdef str charptr_to_str(const char* s):
#     if s == NULL:
#         return None
#     return s.decode(TEXT_ENCODING, ERROR_HANDLER)