from .FileIO cimport *

cdef object getAssembleArray(dict chromCltData)
cpdef wtdbg2Assemble(int[:, :] assembleArray, int start, int taskSize, int numThread)