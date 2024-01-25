#ifndef ANNO_UTILS_H
#define ANNO_UTILS_H

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

#endif // ANNO_UTILS_H