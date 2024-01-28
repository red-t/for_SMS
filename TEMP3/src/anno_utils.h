#ifndef ANNO_UTILS_H
#define ANNO_UTILS_H

#include "htslib/faidx.h"
#include "io_utils.h"

/******************
 *** Data Types ***
 ******************/

typedef struct Anno
{
    int     idx;
    int     queryStart;
    int     queryEnd;
    uint8_t strand;
    int     tid;
    int     refStart;
    int     refEnd;
} __attribute__((packed)) Anno;

/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/

/// @brief get all annotate records by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx);

/// @brief init single annotate record
void initAnno(bam1_t *bamRecord, Anno *anno, int idx);

/// @brief annotate polyA/polyT
int annoPolyA(Cluster *cluster, int idx, Anno *annoArray, int numAnno, int leftMost, int rightMost, int leftIdx, int rightIdx);

/// @brief find polyA/polyT region and init single annotate records
int getPolyA(char *flankSeq, int seqLen, int idx, int isA, Anno *annoArray, int numAnno, int rightMost, int leftIdx, int rightIdx);

/// @brief annotate TSD and refine breakpoint for cluster for single cluster
void annoTsd(Cluster *cluster);

/// @brief set TSD and refine breakpoint for cluster
void setTsd(Cluster *cluster, int localStart, int leftEnd, int rightStart);

#endif // ANNO_UTILS_H