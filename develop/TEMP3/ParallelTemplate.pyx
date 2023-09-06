import numpy as np
from .Cluster import build_cluster
from concurrent.futures import ProcessPoolExecutor, as_completed


cdef background_info(str fpath,
                     int nthreads,
                     int *nchroms,
                     float *back_div,
                     float *back_de,
                     float *back_depth,
                     float *back_readlen):
    cdef:
        BamFile rbf = BamFile(fpath, nthreads, "rb")
        int i, tid = 0, maxlen = 0
    
    # Get the longest chromosome
    for i in range(rbf.hdr.n_targets):
        if maxlen < sam_hdr_tid2len(rbf.hdr, i):
            tid = i
            maxlen = sam_hdr_tid2len(rbf.hdr, i)
    
    cdef:
        Iterator ite = Iterator(rbf, tid)
        int retval, N = 0, M = 499980, alnlen
        object divs, des, readlens, templatef, templatei
        float[::1] divs_view
        float[::1] des_view
        int[::1] readlens_view
        int64_t sum_alnlen = 0
        float div

    templatef   = np.zeros(500000, dtype=np.float32)
    templatei   = np.zeros(500000, dtype=np.int32)
    divs        = np.zeros(500000, dtype=np.float32)
    des         = np.zeros(500000, dtype=np.float32)
    readlens    = np.zeros(500000, dtype=np.int32)
    divs_view       = divs
    des_view        = des
    readlens_view   = readlens

    # Estiamte background divergence & coverage & readlen
    while 1:
        if N > M:
            divs        = np.concatenate((divs, templatef))
            des         = np.concatenate((des, templatef))
            readlens    = np.concatenate((readlens, templatei))
            divs_view       = divs
            des_view        = des
            readlens_view   = readlens
            M = des.shape[0] - 20

        retval = ite.cnext1()
        if retval > 0:
            if bam_filtered(ite.b):
                continue
            
            get_div(&alnlen, &div, ite.b)
            sum_alnlen      += alnlen
            divs_view[N]     = div
            des_view[N]      = get_de(ite.b)
            readlens_view[N] = ite.b.core.l_qseq
            N += 1
            continue

        nchroms[0]  = rbf.hdr.n_targets
        back_div[0] = np.mean(divs[:N])
        back_de[0]  = np.mean(des[:N])
        back_depth[0]   = float(sum_alnlen) / maxlen
        back_readlen[0] = np.median(readlens[:N])
        del ite; rbf.close(); del rbf
        return
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster_parallel(str fpath,
                                  str rep_path,
                                  str gap_path,
                                  str teref,
                                  str preset,
                                  int nprocess,
                                  int nthreads,
                                  int minl,
                                  int maxdist,
                                  int reftid):
    '''call  in parallel

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
        int nchroms
        float back_div, back_de, back_depth, back_readlen

    # get background information
    background_info(fpath, nprocess, &nchroms, &back_div, &back_de, &back_depth, &back_readlen)
    print("background div: {}\nbackground de: {}\nbackground coverage: {}\nbackground readlen: {}".format(back_div, back_de, back_depth, back_readlen))
    with ProcessPoolExecutor(max_workers=nprocess) as executor:
        futures = set([executor.submit(build_cluster,
                                       fpath,
                                       rep_path,
                                       gap_path,
                                       teref,
                                       preset,
                                       nthreads,
                                       tid,
                                       minl,
                                       maxdist,
                                       back_div,
                                       back_de,
                                       back_depth,
                                       back_readlen) for tid in range(nchroms)])
        # merge result from each process
        for future in as_completed(futures):
            ret = future.result()
            result = {**result, **ret}
    
    return result