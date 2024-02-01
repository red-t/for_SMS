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
void setTrimRegion(Segment *segment, int *start, int *end, int flankSize);


/*************************
 *** Flank Sequence IO ***
 *************************/

/// @brief Data container for reference flank region
typedef struct FlankRegion
{
    int start1;
    int start2;
    int end1;
    int end2;
} FlankRegion;

/// @brief output flank-seq and local-seq for all clusters
void extractRefFlanks(char *refFn, Cluster *cltArray, int startIdx, int endIdx);

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *cluster, FlankRegion *region);

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *cluster, faidx_t *refFa, FlankRegion region);

/// @brief output +-500bp local-seq around cluster position for tsd annotation
void outputLocal(Cluster *cluster, faidx_t *refFa, FlankRegion region);


/*****************************
 *** Insertion Sequence IO ***
 *****************************/
#define bamIsSup(bam) (((bam)->core.flag & BAM_FSUPPLEMENTARY) != 0)
#define isLeftFlank(bam) (strcmp(bam_get_qname((bam)), "0") == 0)
#define isClipInFlank(cigar) (bam_cigar_op((cigar)) == BAM_CSOFT_CLIP && bam_cigar_oplen((cigar)) > 400)
#define isFlankMapped(flag) (((flag) & 120) != 0)
#define isBothFlankMapped(flag) (((flag) & 96) != 0)

/// @brief Data container for insertion region on assembled contig
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

/// @brief output insertion-seq and tsd-containing-seq from contig
void extractIns(Cluster *cluster);

/// @brief define insertion-seq region by Flank-To-Assm alignments
void setInsRegion(int cltTid, int cltIdx, InsRegion *region);

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region);

/// @brief output insertion-seq in FASTA format
void outputInsSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region);

/// @brief get insertion-seq
char *getInsSeq(faidx_t *assmFa, InsRegion region);

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *cluster, faidx_t *assmFa, InsRegion region);

#endif // IO_UTILS_H