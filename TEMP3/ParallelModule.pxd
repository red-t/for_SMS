from .FileIO cimport *
from .Cluster cimport bamIsInvalid, getMapLenAndDiv, getHighQualClts

cpdef object runInParallel(object cmdArgs)