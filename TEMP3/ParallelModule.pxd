from .FileIO cimport *
from .Cluster cimport bamIsInvalid, getMapLenAndDiv
from .Assemble cimport getAssembleArray

cpdef dict runInParallel(object cmdArgs)