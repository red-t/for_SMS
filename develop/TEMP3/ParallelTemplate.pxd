from .AlignmentFileIO cimport BamFile

cpdef dict build_cluster_parallel(str fpath,
                                  str rep_path,
                                  str gap_path,
                                  str teref,
                                  str preset,
                                  int nprocess,
                                  int nthreads,
                                  int minl,
                                  int maxdist)