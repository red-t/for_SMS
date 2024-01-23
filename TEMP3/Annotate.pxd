from .FileIO cimport *
cpdef annotateCluster(Cluster[::1] cltArray, int startIdx, int taskSize, object cmdArgs)