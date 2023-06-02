from .AlignmentFileIO cimport BamFile

cpdef dict build_cluster_parallel(str fpath,
                                  int nprocess,
                                  int nthreads,
                                  int minl,
                                  int maxdist)