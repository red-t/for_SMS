#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H
#include <stdint.h>
#include <math.h>
#include "seg_utils.h"
// #include "AIList.h"
//-------------------------------------------------------------------------------------
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
//-------------------------------------------------------------------------------------
/***********************
 *** Cluster records ***
 ***********************/

/*! @typedef
 @abstract Structure for segment extracted from CIGAR.
 @field  st         0-based start on reference, included
 @field  ed         0-based end on reference, excluded
 @field  st_idx     0-based start index on array of segments, included
 @field  ed_idx     0-based end index on array of segments, excluded, TO BE REMOVED
 @field  nseg       number of segments in this cluster
 @field  strand     the orientation of this insertion
 @field  cloc_flag  bitwise flag of the cluster location
                        1: at normal region
                        2: at repeat/gap boundary
                        4: inside repeat/gap
 @field  ntype      number of segment types in this cluster
 @field  nL         number of left-clip segments in this cluster, TO BE REMOVED
 @field  nM         number of mid-insert segments in this cluster, TO BE REMOVED
 @field  nR         number of right-clip segments in this cluster, TO BE REMOVED
 @field  entropy    entropy computed based on fraction of different segment types
 @field  bratio     balance ratio computed based on number of lef/right-clip segments
 @field  sovh_frac  fraction of segments with short overhang (<100bp)
 @field  lmq_frac   fraction of segments with low mapq (<5)
 @field  dclip_frac fraction of segments with short overhang (<100bp)
 @field  aln1_frac  fraction of segments with loc_flag=1
 @field  aln2_frac  fraction of segments with loc_flag=2
 @field  aln4_frac  fraction of segments with loc_flag=4
 @field  aln8_frac  fraction of segments with loc_flag=8
 @field  aln16_frac fraction of segments with loc_flag=16
 @field  avg_mapq   average mapq of this cluster
 */
typedef struct {
    int32_t     st;
    int32_t     ed;
    int32_t     st_idx;
    int32_t     ed_idx;
    int16_t     nseg;
    uint8_t     strand;
    uint8_t     cloc_flag;
    uint8_t     ntype;
    int16_t     nL;
    int16_t     nM;
    int16_t     nR;
    float       entropy;
    float       bratio;
    float       sovh_frac;
    float       lmq_frac;
    float       dclip_frac;
    float       aln1_frac;
    float       aln2_frac;
    float       aln4_frac;
    float       aln8_frac;
    float       aln16_frac;
    float       avg_mapq;
} __attribute__((packed)) cluster_dtype_struct;


/// Compute number of segment types.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_ntype(cluster_dtype_struct clts[]) {
    clts[0].ntype = (clts[0].ntype&1) + ((clts[0].ntype&2)>>1) + ((clts[0].ntype&4)>>2);
}


/// Compute entropy based on segment types.
/*!
 * @param clts  address to the array of clusters
 */
void clt_entropy(cluster_dtype_struct clts[]);


/// Compute balance ratio based on number of lef/right-clip segments.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_bratio(cluster_dtype_struct clts[]) {
    clts[0].bratio = (MIN(clts[0].nL, clts[0].nR) + 0.01) / (MAX(clts[0].nL, clts[0].nR) + 0.01);
}


/// Compute short overhang fraction.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_sovh(cluster_dtype_struct clts[]) {
    clts[0].sovh_frac = clts[0].sovh_frac / clts[0].nseg;
}


/// Compute low mapq fraction.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_lmq(cluster_dtype_struct clts[]) {
    clts[0].lmq_frac = clts[0].lmq_frac / clts[0].nseg;
}


/// Compute dual-clip alignment fraction.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_dclip(cluster_dtype_struct clts[]) {
    clts[0].dclip_frac = clts[0].dclip_frac / clts[0].nseg;
}


/// Compute fraction of alignment with different loc_flag.
/*!
 * @param clts  address to the array of clusters
 */
void clt_dffloc(cluster_dtype_struct clts[]);


/// Compute average mapq.
/*!
 * @param clts  address to the array of clusters
 */
inline void clt_avgmapq(cluster_dtype_struct clts[]) {
    clts[0].avg_mapq = clts[0].avg_mapq / clts[0].nseg;
}


// /// Compute cluster location flag
// /*!
//  * @param rep_ail	AIList of repeats
//  * @param gap_ail	AIList of gaps
//  * @param clts		address to the array of clusters
//  */
// void clt_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, cluster_dtype_struct clts[]);


// void clt_feat(cluster_dtype_struct clts[], seg_dtype_struct segs[], ailist_t *rep_ail, ailist_t *gap_ail);

#endif // CLUSTER_UTILS_H