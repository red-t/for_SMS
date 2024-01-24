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
void getTrimRegion(Segment *segment, int *startPtr, int *endPtr, int flankSize);

/*************************
 *** Flank Sequence IO ***
 *************************/

/// @brief extract and output flank sequence for all clusters
int outputRefFlankSeqs(char *refFn, Cluster *cltArray, int startIdx, int endIdx);

/*****************************
 *** Insertion Sequence IO ***
 *****************************/
#define bamIsSup(bamRecord) (((bamRecord)->core.flag & BAM_FSUPPLEMENTARY) != 0)
#define isLeftFlank(bamRecord) (strcmp(bam_get_qname((bamRecord)), "0") == 0)

/// @brief extract and output insertion sequence for single cluster
int outputInsSeq(Cluster *cluster);

#endif // IO_UTILS_H