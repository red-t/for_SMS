#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <string.h>
#include "htslib/faidx.h"
#include "cluster_utils.h"

int outputRefFlankSeqs(char *refFn, Cluster *cltArray, int startIdx, int endIdx);

#endif // IO_UTILS_H