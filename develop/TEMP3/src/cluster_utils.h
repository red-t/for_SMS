#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H
#include <stdint.h>
#include <math.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "seg_utils.h"
#include "AIList.h"
//-------------------------------------------------------------------------------------
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
//-------------------------------------------------------------------------------------
/***********************
 *** Cluster records ***
 ***********************/

/*! @typedef
 @abstract Structure cluster merged from segments.
 @field  st         0-based start on reference, included
 @field  ed         0-based end on reference, excluded
 @field  st_idx     0-based start index on array of segments, included
 @field  ed_idx     0-based end index on array of segments, excluded, TO BE REMOVED
 @field  nseg       number of segments in this cluster (normalized by estimated background coverage)
 @field  strand     the orientation of this insertion
                        1: forward
                        2: reverse
 @field  single     whether the cluster only have single support read
                        0: multiple support reads
                        1: single read & 1 alignment
                        2: single read & 2 alignments
 @field  cloc_flag  bitwise flag of the cluster location
                        1: at normal region
                        2: at repeat/gap boundary
                        4: inside repeat/gap
 @field  ntype      number of segment types in this cluster
 @field  entropy    entropy computed based on fraction of different segment types
 @field  bratio     balance ratio computed based on number of lef/right-clip segments
 @field  lmq_frac   fraction of segments with low mapq (<5)
 @field  dclip_frac fraction of "dual-clip" alignments
 @field  aln1_frac  fraction of segments with loc_flag=1
 @field  aln2_frac  fraction of segments with loc_flag=2
 @field  aln4_frac  fraction of segments with loc_flag=4
 @field  aln8_frac  fraction of segments with loc_flag=8
 @field  aln16_frac fraction of segments with loc_flag=16
 @field  avg_mapq   average mapq of this cluster
 @field  avg_AS     average per base alignment score
 @field  avg_qfrac  average query aligned fraction
 @field  avg_div    average per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
 @field  avg_de     average per-base gap-compressed divergence (normalized background)
 @field  back_div   background div
 @field  back_de    background de
 @field  back_depth background depth
 @field  back_readlen   background read length
 */
typedef struct {
    int32_t     st;
    int32_t     ed;
    int32_t     st_idx;
    int32_t     ed_idx;
    float_t     nseg;
    uint16_t    strand;
    uint8_t     single;
    uint8_t     cloc_flag;
    uint8_t     ntype;
    float_t     entropy;
    float_t     bratio;
    float_t     lmq_frac;
    float_t     dclip_frac;
    float_t     aln1_frac;
    float_t     aln2_frac;
    float_t     aln4_frac;
    float_t     aln8_frac;
    float_t     aln16_frac;
    float_t     avg_mapq;
    float_t     avg_AS;
    float_t     avg_qfrac;
    float_t     avg_div;
    float_t     avg_de;
    float_t     back_div;
    float_t     back_de;
    float_t     back_depth;
    float_t     back_readlen;
} __attribute__((packed)) cluster_dtype_struct;


/// Compute cluster location flag
/*!
 * @param rep_ail	AIList of repeats
 * @param gap_ail	AIList of gaps
 * @param clts		address to the cluster record
 */
void clt_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, cluster_dtype_struct clts[]);


/// Compute entropy based on segment types.
/*!
 * @param nL    number of left-clip segments
 * @param nM    number of mid-insert segments
 * @param nR    number of right-clip segments
 * @param nseg  number of total segments in the cluster
 */
float_t clt_entropy(int32_t nL, int32_t nM, int32_t nR, int32_t nseg);


/// Compute fraction of alignment with different loc_flag.
/*!
 * @param clts  address to the cluster record
 * @param nseg  number of segments in the cluster
 */
void clt_dffloc(cluster_dtype_struct clts[], int32_t nseg);


/// Compute features of a cluster record.
/*!
 * @param clts      address to the cluster record
 * @param segs      address to the arrary of segments
 * @param repail    AIList of repeats
 * @param gapail    AIList of gaps
 * @param minovh    minimum length of segment overhang
 */
void cclt_feat(cluster_dtype_struct clts[],
               seg_dtype_struct segs[],
               ailist_t *rep_ail,
               ailist_t *gap_ail,
               float_t back_div,
               float_t back_de,
               float_t back_depth,
               float_t back_readlen,
               int minovh,
               int tid,
               htsFile *htsfp,
               bam1_t *b1,
               bam1_t *b2);

#endif // CLUSTER_UTILS_H