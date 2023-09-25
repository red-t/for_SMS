#ifndef SEG_UTILS_H
#define SEG_UTILS_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "AIList.h"

/**********************
 *** CIGAR resolver ***
 **********************/

/// Get whether the CIGAR operation is clip/insert.
/*!
 * @param op   CIGAR operation
 * @return     1 if the CIGAR operation is clip/insert, 0 if not
 */
inline int is_clip_or_insert(uint32_t op) {
    return op == BAM_CSOFT_CLIP || op == BAM_CINS;
}


/// Get whether the CIGAR operation is match/equal/diff.
/*!
 * @param op   CIGAR operation
 * @return     1 if the CIGAR operation is match/equal/diff, 0 if not
 */
inline int is_match(uint32_t op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}


/// Get whether the CIGAR operation is del/skip.
/*!
 * @param op   CIGAR operation
 * @return     1 if the CIGAR operation is del/skip, 0 if not
 */
inline int is_del_or_skip(uint32_t op) {
    return op == BAM_CDEL || op == BAM_CREF_SKIP;
}


/***********************
 *** Segment records ***
 ***********************/

#define LEFT_CLIP   1
#define MID_INSERT  2
#define RIGHT_CLIP  4
#define DUAL_CLIP   5

/*! @typedef
 @abstract Structure for segment extracted from CIGAR.
 @field  flag       bitwise flag of the alignment
 @field  mapq       mapping quality
 @field  qst        0-based start on query sequence, included
 @field  qed        0-based end on query sequence, excluded
 @field  pos        0-based position on reference, included
 @field  sflag      bitwise flag of the segment type
                        1: left clip
                        2: internal insert
                        4: right clip
 @field  rflag      bitwise flag of the read/alignment type, that is, combination of sflag(s)
 @field  offset     file offset of the corresponding alignment record
 @field  refst      0-based reference start of the corresponding alignment record
 @field  refed      0-based reference end of the corresponding alignment record
 @field  ith        indicates that this is the i-th segments from the same alignments
 @field  nseg       number of segments extracted from the same alignment
 @field  overhang   number of mathced bases on one side of the segment breakpoint
 @field  nmatch     number of mathced bases of the alignment
 @field  loc_flag   bitwise flag of the alignment location
                        1:  at least one side at normal region
                        2:  both sides inside repeat/gap
                        4:  both sides at repeat/gap boundary
                        8:  one side at repeat/gap boundary, the other at normal region
                        16: one side at repeat/gap boundary, the other inside repeat/gap
 @field  nmap       number of alignments mapped to TE
 @field  lmap       summary query length of the TE alignments
 @field  sumAS      summary per base alignment score of the TE alignments
 @field  sumdiv     summary per-base divergence of the TE alignments
 @field  sumde      summary gap-compressed per-base divergence("de") of the TE alignments
 @field  cnst       number of mathced bases of the alignment
 @field  st_idx     0-based start index on array of TE alignments, included
 @field  ed_idx     0-based end index on array of TE alignments, excluded
 @field  TE         majority TE-tid of the segment
 */
typedef struct {
    uint16_t    flag;
    uint8_t     mapq;
    int32_t     qst;
    int32_t     qed;
    int32_t     rpos;
    uint8_t     sflag;
    uint8_t     rflag;
    int64_t     offset;
    int32_t     refst;
    int32_t     refed;
    uint8_t     ith;
    uint8_t     nseg;
    int32_t     overhang;
    int32_t     nmatch;
    uint8_t     loc_flag;
    uint8_t     nmap;
    int32_t     lmap;
    float_t     sumAS;
    float_t     sumdiv;
    float_t     sumde;
    uint16_t    cnst;
    int32_t     st_idx;
    int32_t     ed_idx;
    int32_t     TE;
} __attribute__((packed)) seg_dtype_struct;


/*! @typedef
 @abstract Structure for TE alignment.
 @field  idx    index of the source segment record
 @field  AS     alignment score
 @field  qst    query start (original orientation)
 @field  qed    query end (original orientation)
 @field  alnlen alignment length
 @field  div    per-base divergence
 @field  de     gap-compressed per-base divergence ("de")
 @field  flag   bitwise flag of the alignment
 @field  TE     tid of the TE alignment
 */
typedef struct {
    int32_t idx;
    int32_t AS;
    int32_t qst;
    int32_t qed;
    int32_t alnlen;
    float_t div;
    float_t de;
    int16_t flag;
    int32_t TE;
} __attribute__((packed)) tealn_dtype_struct;


/// parse alignment's CIGAR and extract segments.
/*!
 * @param bam		Alignment record
 * @param segs      address to the segments arrary record
 * @param offset    offset of the alignments in BAM file
 * @param minl      minimum length of segment
 */
int parse_cigar(bam1_t *bam, seg_dtype_struct segs[], int64_t offset, int minl);


/// Get whether the alignment is dual-clip.
/*!
 * @param rflag     rflag of the segment, representing type of the corresponding alignment
 * @return          1 if the alignment is dual-clip, 0 if not
 */
inline int aln_is_dclip(uint8_t rflag) {
    return (rflag & DUAL_CLIP) == 5;
}


/// Get whether the query is secondary alignment by flag.
/*!
 * @param flag  bitwise flag of the query alignment
 * @return      1 if query is secondary, 0 if not
 */
inline int aln_is_second(uint16_t flag) {
    return (flag & BAM_FSECONDARY) != 0;
}


/// Get whether the query is on the reverse strand by flag.
/*!
 * @param flag  bitwise flag of the query alignment
 * @return      1 if query is on the reverse strand, 0 if not
 */
inline int aln_is_rev(uint16_t flag) {
    return (flag & BAM_FREVERSE) != 0;
}


/// Update overhang for mid-insert type segment
/*!
 * @param overhang  original overhang of the segment
 * @param nmatch    nummber of match bases of the corresponding alignment
 * @return          Updated overhang, <= original overhang
 */
inline int32_t update_overhang(int32_t overhang, int32_t nmatch) {
    return (nmatch - overhang) < overhang ? (nmatch - overhang) : overhang;
}


/// Compute alignment location flag
/*!
 * @param rep_ail	AIList of repeats
 * @param gap_ail	AIList of gaps
 * @param segs		address to the segment record
 */
void aln_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, seg_dtype_struct segs[]);


/// Compute features of a segment record
/*!
 * @param segs		address to the segment record
 * @param rep_ail	AIList of repeats
 * @param gap_ail	AIList of gaps
 */
void cseg_feat(seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail);


/// Compute features of a segment record from TE alignment
/*!
 * @param segs		address to the segment record
 * @param tealns    address to the TE alignments array
 * @param i         index of the used alignment record
 */
void cseg_feat_te(seg_dtype_struct segs[], tealn_dtype_struct tealns[], int i);


/// Determine the majority TE type of a segment record
/*!
 * @param segs		address to the segment record
 * @param tealns    address to the TE alignments array
 * @param TEs       address to the TE frequency array
 * @param TE_size   size of the TE frequency array
 */
void cseg_tetype(seg_dtype_struct segs[], tealn_dtype_struct tealns[], int TEs[], int TE_size);


/*************************
 *** Alignment records ***
 *************************/

/// Get whether the query is secondary alignment
/*!
 @param  b  pointer to an alignment
 @return    1 if query is secondary, 0 if not
 */
#define bam_is_second(b) (((b)->core.flag & BAM_FSECONDARY) != 0)

/// Get whether the query is secondary or unmapped
/*!
 @param  b  pointer to an alignment
 @return    1 if query is secondary or unmapped, 0 if not
 */
#define bam_filtered(b) (((b)->core.flag & BAM_FSECONDARY) != 0 || ((b)->core.flag & BAM_FUNMAP) != 0)

/// Get "per-base divergence" of the alignment
/*!
 @param  b  pointer to an alignment
 @return    per-base divergence
 */
void get_div(int32_t *alnlen, float_t *div, bam1_t *b);

/// Get "gap-compressed per-base divergence" of the alignment
/*!
 @param  b  pointer to an alignment
 @return    gap-compressed per-base divergence
 */
#define get_de(b) (bam_aux2f(bam_aux_get((b), "de")))

/// Set alignment record memory policy
/**
   @param b       Alignment record
   @param policy  Desired policy
*/
static inline void bam_set_mempolicy1(bam1_t *b, uint32_t policy) {
    b->mempolicy = policy;
}

/// Get alignment record memory policy
/** @param b    Alignment record

    See bam_set_mempolicy()
 */
static inline uint32_t bam_get_mempolicy1(bam1_t *b) {
    return b->mempolicy;
}

/// bam1_t data (re)allocation
int sam_realloc_bam_data1(bam1_t *b, size_t desired);

static inline int realloc_bam_data1(bam1_t *b, size_t desired)
{
    if (desired <= b->m_data) return 0;
    return sam_realloc_bam_data1(b, desired);
}

/// copy sequence from source alignment to destination alignment
/**
   @param src   source alignment record     
   @param dest  destination alignment record
   @param idx   index of the segment in the arrary, will be used as qname of dest
   @param qst   query start, will copy sequence from src from qst
   @param qed   query end
   @return      dest data length if success, -1 if failed
*/
int bam_trim1(bam1_t *src, bam1_t *dest, int32_t idx, int32_t qst, int32_t qed);


/// parse transposon alignment and record record information with array
/*!
 * @param bam		Alignment record
 * @param tealns	address to the te alignments arrary record
 */
void parse_tealns(bam1_t *bam, tealn_dtype_struct tealns[]);

#endif // SEG_UTILS_H