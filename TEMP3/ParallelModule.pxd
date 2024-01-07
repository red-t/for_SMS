from .FileIO cimport BamFile, Iterator
from .Cluster cimport bamIsInvalid, getMapLenAndDiv
from .HtslibExternal cimport *

cpdef dict runInParallel(object cmdArgs)