from .FileIO cimport *
from .Cluster cimport bamIsInvalid, getMapLenAndDiv, getAnnoArray
from .Assemble cimport getAssembleArray

cpdef dict runInParallel(object cmdArgs)