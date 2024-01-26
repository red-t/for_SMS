#ifndef ANNO_UTILS_H
#define ANNO_UTILS_H

#include "htslib/faidx.h"
#include "cluster_utils.h"

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

#endif // ANNO_UTILS_H