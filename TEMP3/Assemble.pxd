from .FileIO cimport *
from cpython cimport PyBytes_FromStringAndSize

cpdef assembleCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs)