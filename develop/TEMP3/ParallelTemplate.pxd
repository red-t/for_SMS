from .AlignmentFileIO cimport BamFile, Iterator
from .Cluster cimport bamIsInvalid, getMapLenAndDiv, getDe
from .htslib_external cimport *

cpdef dict buildClusterParallel(object args)