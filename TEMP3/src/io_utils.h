#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <string.h>
#include "htslib/faidx.h"
#include "cluster_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/
#define isLowQualClt(cluster) (((cluster)->isInBlacklist) || ((cluster)->probability <= 0.5))
#define isSomaClt(cluster) ((cluster)->cltType != 0)

/// @brief select a segment to output from cluster
int getOuputSegIdx(Cluster *cluster, Segment *segArray, Args args);

/// @brief get extended region of the segment
void setTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize);

/*************************
 *** Flank Sequence IO ***
 *************************/

typedef struct FlankRegion
{
    int start1;
    int start2;
    int end1;
    int end2;
} FlankRegion;

/// @brief extract and output flank sequence for all clusters
void extractRefFlanks(char *refFn, Cluster *cltArray, int startIdx, int endIdx);

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *cluster, FlankRegion *region);

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *cluster, faidx_t *refFa, FlankRegion region);

/// @brief output +-500bp around cluster position for tsd annotation
void outputRefLocal(Cluster *cluster, faidx_t *refFa, FlankRegion region);

/*****************************
 *** Insertion Sequence IO ***
 *****************************/
#define bamIsSup(bamRecord) (((bamRecord)->core.flag & BAM_FSUPPLEMENTARY) != 0)
#define isLeftFlank(bamRecord) (strcmp(bam_get_qname((bamRecord)), "0") == 0)
#define isClipInFlank(cigar) (bam_cigar_op((cigar)) == BAM_CSOFT_CLIP && bam_cigar_oplen((cigar)) > 400)
#define isFlankMapped(flag) (((flag) & 120) != 0)
#define isBothFlankMapped(flag) (((flag) & 96) != 0)

    typedef struct InsRegion
{
    int start2;
    int end1;
    int tid1;
    int tid2;
    int len1;
    uint32_t cigar1;
    uint32_t cigar2;
    uint16_t flag;
} InsRegion;

/// @brief extract and output insertion-seq and flank-seq
void extractIns(Cluster *cluster);

/// @brief define insertion sequence region
void setInsRegion(int cltTid, int cltIdx, InsRegion *region);

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region);

/// @brief output insertion sequence in FASTA format
void outputInsSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region);

/// @brief get insertion sequence
char *getInsSeq(faidx_t *assmFa, InsRegion region);

/// @brief output flank sequences from assembly for tsd annotation
void outputTsdSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region);

#endif // IO_UTILS_H