#ifndef ANNO_UTILS_H
#define ANNO_UTILS_H

#include "htslib/faidx.h"
#include "io_utils.h"

/******************
 *** Structures ***
 ******************/

/// @brief Data container for annotation record
typedef struct Anno
{
    int     idx;
    int     cltTid;
    int     cltIdx;
    int     queryStart;
    int     queryEnd;
    uint8_t strand;
    int     tid;
    int     refStart;
    int     refEnd;
} __attribute__((packed)) Anno;

/// @brief Data container for polyA
typedef struct PolyA
{
    int idx;
    int cltTid;
    int cltIdx;
    int leftMost;
    int leftIdx;
    int rightMost;
    int rightIdx;
    int isA;
    int seqLen;
} PolyA;


/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *cluster, Anno *annoArray, int idx);

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *cluster, Anno *annoArray, int numAnno, PolyA polyA);

/// @brief Find and record single polyA/polyT
int getPolyA(char *flankSeq, Anno *annoArray, int numAnno, PolyA polyA);

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annoTsd(Cluster *cluster);

/// @brief Find TSD and refine breakpoint
void setTsd(Cluster *cluster, int localStart, int leftEnd, int rightStart);

/// @brief Output annotation records
void outPutAnno(Anno *annoArray, int numAnno, const char *outFn);

#endif // ANNO_UTILS_H