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

    # Compute features of a cluster record
    # @param clts		address to the cluster record
    # @param segs       address to the arrary of segments
    # @param rep_ail	AIList of repeats
    # @param gap_ail    AIList of gaps
    void clt_feat(cluster_dtype_struct *, seg_dtype_struct *, ailist_t *rep_ail, ailist_t *gap_ail)
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
    # @return      0 if query is secondary, 1 if not
    int aln_not_second(uint16_t flag)

    # Get whether the query is on the reverse strand by flag.
    # @param flag  bitwise flag of the query alignment
    # @return      1 if query is on the reverse strand, 0 if not
    int aln_is_rev(uint16_t flag)

    # Compute features of a segment record
    # @param segs		address to the segment record
    # @param rep_ail	AIList of repeats
    # @param gap_ail	AIList of gaps
    void seg_feat(seg_dtype_struct *, ailist_t *rep_ail, ailist_t *gap_ail)

    #########################
    ### Alignment records ###
    #########################

    # Compute features of a segment record
    # @param src   source alignment record     
    # @param dest  destination alignment record
    # @param idx   index of the segment in the arrary, will be used as qname of dest
    # @param qst   query start, will copy sequence from src from qst
    # @param qed   query end
    # @return      dest data length if success, -1 if failed
    int bam_trim1(bam1_t *src, bam1_t *dest, int32_t idx, int32_t qst, int32_t qed)
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