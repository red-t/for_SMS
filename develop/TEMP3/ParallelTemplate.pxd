from .AlignmentFileIO cimport BamFile, Iterator
from .Cluster cimport bam_filtered, get_div, get_de
from .htslib_external cimport *

cpdef dict build_cluster_parallel(str fpath,
                                  str rep_path,
                                  str gap_path,
                                  str teref,
                                  str preset,
                                  int nprocess,
                                  int nthreads,
                                  int minl,
                                  int maxdist,
                                  int reftid)