#ifndef TEMP_UTIL_H
#define TEMP_UTIL_H

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

#define LEFT_CLIP       0x1
#define MID_INSERT      0x2
#define RIGHT_CLIP      0x4
#define DUAL_CLIP       0x5

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


#endif // TEMP_UTIL_H