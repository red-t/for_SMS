from .AlignmentFileIO cimport BamFile, Iterator
from .SegmentParser cimport seg_dtype_struct
from .htslib_external cimport *

cdef extern from "src/cluster_utils.h" nogil:
    #######################
    ### Cluster records ###
    #######################
    ctypedef packed struct cluster_dtype_struct:
        int32_t     st
        int32_t     ed
        int32_t     st_idx
        int32_t     ed_idx
        int16_t     nseg
        uint8_t     strand
        uint8_t     cloc_flag
        uint8_t     ntype
        int16_t     nL
        int16_t     nM
        int16_t     nR
        float       entropy
        float       bratio
        float       sovh_frac
        float       lmq_frac
        float       dclip_frac
        float       aln1_frac
        float       aln2_frac
        float       aln4_frac
        float       aln8_frac
        float       aln16_frac
        float       avg_mapq
    
    # Compute number of segment types.
    # @param clts  address to the array of clusters
    void clt_ntype(cluster_dtype_struct *)

    # Compute entropy based on segment types.
    # @param clts  address to the array of clusters
    void clt_entropy(cluster_dtype_struct *)

    # Compute balance ratio based on segment types.
    # @param clts  address to the array of clusters
    void clt_bratio(cluster_dtype_struct *)

    # Compute short overhang fraction.
    # @param clts  address to the array of clusters
    void clt_sovh(cluster_dtype_struct *)

    # Compute low mapq fraction.
    # @param clts  address to the array of clusters
    void clt_lmq(cluster_dtype_struct *)

    # Compute dual-clip alignment fraction.
    # @param clts  address to the array of clusters
    void clt_dclip(cluster_dtype_struct *)

    # Compute fraction of alignment with different loc_flag.
    # @param clts  address to the array of clusters
    void clt_dffloc(cluster_dtype_struct *)

    # Compute average mapq.
    # @param clts  address to the array of clusters
    void clt_avgmapq(cluster_dtype_struct *)

    # void clt_feat(cluster_dtype_struct *, seg_dtype_struct *, ailist_t *rep_ail, ailist_t *gap_ail)
#
# ---------------------------------------------------------------
#
cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################

    uint8_t LEFT_CLIP
    uint8_t MID_INSERT
    uint8_t RIGHT_CLIP
    uint8_t DUAL_CLIP

    # Get whether the alignment is dual-clip.
    # @param rflag  rflag of the segment, representing type of the corresponding alignment
    # @return       1 if the alignment is dual-clip, 0 if not
    int aln_is_dclip(uint8_t rflag)

    # Get whether the query is secondary alignment by flag.
    # @param flag  bitwise flag of the query alignment
    # @return      1 if query is secondary, 0 if not
    int aln_is_second(uint16_t flag)

    # Get whether the query is on the reverse strand by flag.
    # @param flag  bitwise flag of the query alignment
    # @return      1 if query is on the reverse strand, 0 if not
    int aln_is_rev(uint16_t flag)

    # Update overhang for mid-insert type segment
    # @param overhang  original overhang of the segment
    # @param nmatch    nummber of match bases of the corresponding alignment
    # @return          Updated overhang, <= original overhang
    int32_t update_overhang(int32_t overhang, int32_t nmatch)
#
# ---------------------------------------------------------------
#
cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    ctypedef struct ctg_t:
        pass

    ctypedef struct ailist_t:
        ctg_t *ctg

    # Initialize ailist_t
    ailist_t *ailist_init()

    # Free ailist data
    void ailist_destroy(ailist_t *ail)

    # Add intervals from BED-like file
    # @param fn   filename of the BED-like file
    # @param chr  name of the specified chromosome, only intervals
    #             on this chromosome will be added
    void readBED(ailist_t *ail, const char* fn, char* chr)

    # Construct ailist: decomposition and augmentation
    # @param cLen minimum coverage length, default=20
    void ailist_construct(ailist_t *ail, int cLen)

    # Compute alignment location flag
    # @param rep_ail	AIList of repeats
    # @param gap_ail	AIList of gaps
    # @param segs		address to the array of segments
    void aln_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, seg_dtype_struct *)

    # Compute cluster location flag
    # @param rep_ail	AIList of repeats
    # @param gap_ail	AIList of gaps
    # @param clts		address to the array of clusters
    void clt_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, cluster_dtype_struct *)

#
# ---------------------------------------------------------------
#
cdef object merge_segments(seg_dtype_struct[::1] segs,
                                            int maxdist)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str fpath,
                         str rep_path,
                         str gap_path,
                         int threads,
                         int tid,
                         int minl,
                         int maxdist)