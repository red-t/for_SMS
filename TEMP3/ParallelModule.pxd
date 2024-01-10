from .FileIO cimport *
from .Cluster cimport bamIsInvalid, getMapLenAndDiv, getHighQualClts

cpdef dict runInParallel(object cmdArgs)