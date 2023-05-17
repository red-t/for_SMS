# filename: ParallelTemplate.pyx

from pysam import AlignmentFile
from Cluster import build_cluster
from concurrent.futures import ProcessPoolExecutor, as_completed


cpdef dict build_cluster_parallel(str fpath,
                                  int nprocess,
                                  int nthreads,
                                  int minl,
                                  int maxdist):
    '''call build_cluster in parallel

    call build_cluster in parallel with multiple process.

    Parameters:
    -----------
        fpath: str
            path of input BAM file.

        nprocess: int
            maximum nummber of worker processes.
        
        nthreads: int
            maximum nummber of extra threads to use in each sub 
            process.
        
        minl: int
            minimun segment length, cigar operation with length < 
            minl will not be used to create insert segment.
        
        maxdist: int
            max merging distance, segments with distance larger t-
            han maxdist will not be merged in to the same cluster.

    Returns:
    --------
        result: dict
            dictionary of clusters, tid -> list_of_Cluster.
    '''
    cdef set    futures
    cdef list   idx
    cdef dict   result={}, ret
    cdef int    nchroms, i
    
    bf = AlignmentFile(fpath, "rb")
    nchroms = bf.nreferences
    bf.close(); del bf

    with ProcessPoolExecutor(max_workers=nprocess) as executor:
        futures = set([executor.submit(build_cluster, fpath,
                                       nthreads, i,
                                       minl, maxdist) for i in range(nchroms)])
        # merge result from each process
        for future in as_completed(futures):
            ret = future.result()
            result = {**result, **ret}
    
    return result