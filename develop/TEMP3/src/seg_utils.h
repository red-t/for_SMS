#ifndef SEG_UTILS_H
#define SEG_UTILS_H

#include <stdint.h>
#include "htslib/sam.h"

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
 @field  lqseq      length of the query sequence (read)
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
 */
typedef struct {
    uint16_t    flag;
    uint8_t     mapq;
    int32_t     qst;
    int32_t     qed;
    int32_t     rpos;
    int32_t     lqseq;
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
} __attribute__((packed)) seg_dtype_struct;

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


/*************************
 *** Alignment records ***
 *************************/

/// Get whether the query is secondary alignment
/*!
 @param  b  pointer to an alignment
 @return    1 if query is secondary, 0 if not
 */
#define bam_is_second(b) (((b)->core.flag & BAM_FSECONDARY) != 0)


#endif // SEG_UTILS_H