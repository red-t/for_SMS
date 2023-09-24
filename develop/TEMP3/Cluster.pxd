from .AlignmentFileIO cimport BamFile, Iterator
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
        float       nseg
        uint16_t    strand
        uint8_t     single
        uint8_t     cloc_flag
        uint8_t     ntype
        float       entropy
        float       bratio
        float       lmq_frac
        float       dclip_frac
        float       aln1_frac
        float       aln2_frac
        float       aln4_frac
        float       aln8_frac
        float       aln16_frac
        float       avg_mapq
        float       avg_AS
        float       avg_qfrac
        float       avg_div
        float       avg_de
        float       back_div
        float       back_de
        float       back_depth
        float       back_readlen
        int16_t     flag

    # Compute features of a cluster record
    # @param clts           address to the cluster record
    # @param segs           address to the arrary of segments
    # @param rep_ail        AIList of repeats
    # @param gap_ail        AIList of gaps
    # @param back_de        estimated background gap-compressed per-base divergence
    # @param back_depth     estimated background coverage
    # @param back_readlen   estimated background read length
    # @param minovh         minimum length of segment overhang
    # @param tid            target id
    # @param htsfp          pointer to the htsFile
    # @param idx            pointer to the hts_idx_t
    # @param b1             first bam1_t record
    # @param b2             second bam1_t record
    void cclt_feat(cluster_dtype_struct *,
                   seg_dtype_struct *,
                   ailist_t *rep_ail,
                   ailist_t *gap_ail,
                   float back_div,
                   float back_de,
                   float back_depth,
                   float back_readlen,
                   int minovh,
                   int tid,
                   htsFile *htsfp,
                   bam1_t *b1,
                   bam1_t *b2)
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

    ctypedef packed struct seg_dtype_struct:
        uint16_t    flag
        uint8_t     mapq
        int32_t     qst
        int32_t     qed
        int32_t     rpos
        uint8_t     sflag
        uint8_t     rflag
        int64_t     offset
        int32_t     refst
        int32_t     refed
        uint8_t     ith
        uint8_t     nseg
        int32_t     overhang
        int32_t     nmatch
        uint8_t     loc_flag
        uint8_t     nmap
        int32_t     lmap
        float       sumAS
        float       sumdiv
        float       sumde
        uint16_t    cnst

    ctypedef packed struct tealn_dtype_struct:
        int32_t idx
        int32_t AS
        int32_t qst
        int32_t qed
        int32_t alnlen
        float   div
        float   de
        int16_t flag
    
    # parse alignment's CIGAR and extract segments.
    # @param bam    Alignment record
    # @param segs   address to the segments arrary record
    # @param offset offset of the alignments in BAM file
    # @param minl   minimum length of segment
    int parse_cigar(bam1_t *bam, seg_dtype_struct *, int64_t offset, int minl)

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

    # Compute features of a segment record
    # @param segs		address to the segment record
    # @param rep_ail	AIList of repeats
    # @param gap_ail	AIList of gaps
    void cseg_feat(seg_dtype_struct *, ailist_t *rep_ail, ailist_t *gap_ail)

    # Compute features of a segment record from TE alignment
    # @param segs		address to the segment record
    # @param tealns    address to the TE alignments array
    # @param i         index of the used alignment record
    void cseg_feat_te(seg_dtype_struct *, tealn_dtype_struct *, int i)

    #########################
    ### Alignment records ###
    #########################

    # Get whether the query is secondary or unmapped
    # @param  b  pointer to an alignment
    # @return    1 if query is secondary or unmapped, 0 if not
    int bam_filtered(bam1_t *b)

    # Get "per-base divergence" of the alignment
    # @param  b  pointer to an alignment
    # @return    per-base divergence
    void get_div(int *alnlen, float *div, bam1_t *b)
    
    # Get "gap-compressed per-base divergence" of the alignment
    # @param  b  pointer to an alignment
    # @return    gap-compressed per-base divergence
    float get_de(bam1_t *b)

    # Compute features of a segment record
    # @param src   source alignment record     
    # @param dest  destination alignment record
    # @param idx   index of the segment in the arrary, will be used as qname of dest
    # @param qst   query start, will copy sequence from src from qst
    # @param qed   query end
    # @return      dest data length if success, -1 if failed
    int bam_trim1(bam1_t *src, bam1_t *dest, int32_t idx, int32_t qst, int32_t qed)

    # parse transposon alignment and record record information with array
    # @param bam    Alignment record
    # @param tealns address to the arrary record
    void parse_tealns(bam1_t *bam, tealn_dtype_struct *)
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
cdef object merge_seg(seg_dtype_struct[::1] segs,
                      int maxdist)
#
# ---------------------------------------------------------------
#
cpdef dict build_cluster(str fpath,
                         str rep_path,
                         str gap_path,
                         str teref,
                         str preset,
                         int threads,
                         int tid,
                         int minl,
                         int maxdist,
                         float back_div,
                         float back_de,
                         float back_depth,
                         float back_readlen)