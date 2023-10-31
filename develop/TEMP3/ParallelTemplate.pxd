from .AlignmentFileIO cimport BamFile, Iterator
from .Cluster cimport bamIsInvalid, getMapLenAndDiv
from .htslib_external cimport *

cpdef dict buildClusterParallel(object cmdArgs)