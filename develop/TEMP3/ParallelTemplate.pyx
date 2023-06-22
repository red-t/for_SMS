from .Cluster import build_cluster
from concurrent.futures import ProcessPoolExecutor, as_completed


cdef float background_div(int N,
                          Iterator ite):
    cdef:
        int n=0, retval
        float div=0
    
    while n < N:
        retval = ite.cnext2()
        if retval > 0:
            if bam_filtered(ite.b):
                continue
            
            div += get_div(ite.b)
            n += 1
            continue

        div = div/n
        return div
    
    div = div/n
    return div


cdef background_info(str fpath,
                     int nthreads,
                     int *nchroms,
                     float *bdiv,
                     float *coverage):
    cdef:
        BamFile rbf = BamFile(fpath, nthreads, "rb")
        int i, tid=0
        int maxlen = 0
    
    for i in range(rbf.hdr.n_targets):
        if maxlen < sam_hdr_tid2len(rbf.hdr, i):
            tid = i
            maxlen = sam_hdr_tid2len(rbf.hdr, i)

    cdef:
        Iterator ite = Iterator(rbf, tid)
        int64_t sumlen = 0
        int n = 0, retval
        float div = 0
    
    while 1:
        retval = ite.cnext1()
        if retval > 0:
            if bam_filtered(ite.b):
                continue
            
            div += get_div(ite.b)
            sumlen += ite.b.core.l_qseq
            n += 1
            continue

        nchroms[0] = rbf.hdr.n_targets
        bdiv[0] = div/n
        coverage[0] = sumlen/maxlen
        del ite; rbf.close(); del rbf
        return




cpdef dict build_cluster_parallel(str fpath,
                                  str rep_path,
                                  str gap_path,
                                  str teref,
                                  str preset,
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
    cdef:
        dict result={}, ret
        set futures
        int nchroms, i
        float div
        BamFile bf = BamFile(fpath, nthreads, "rb")
        Iterator ite = Iterator(bf)
    
    # get number of target chromosomes
    nchroms = bf.hdr.n_targets
    div = background_div(50000, ite)
    print("background divergence: ", div)
    del ite; bf.close(); del bf

    with ProcessPoolExecutor(max_workers=nprocess) as executor:
        futures = set([executor.submit(build_cluster, fpath, rep_path, gap_path, teref, preset,
                                       nthreads, i, minl, maxdist, div) for i in range(nchroms)])
        # merge result from each process
        for future in as_completed(futures):
            ret = future.result()
            result = {**result, **ret}
    
    return result