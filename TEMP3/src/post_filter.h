#ifndef POST_FILTER_H
#define POST_FILTER_H

#include "anno_utils.h"

// Perform post-filtering for different TE class
void postFilter(Cluster *clt);

// Perform post-filtering for LINE, SINE, RETROPOSON
void filterLINE(Cluster *clt);

// Perform post-filtering for LTR, DNA
void filterLTR(Cluster *clt);

#endif // POST_FILTER_H