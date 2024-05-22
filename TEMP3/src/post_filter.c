#include "post_filter.h"

// Perform post-filtering for different TE class
void postFilter(Cluster *clt)
{
    if ((clt->flag & (CLT_LINE|CLT_SINE|CLT_RETROPOSON)) != 0)
        filterLINE(clt);
    
    if ((clt->flag & CLT_LTR) != 0)
        filterLTR(clt);

    if ((clt->flag & CLT_DNA) != 0)
        filterDNA(clt);
}

// Perform post-filtering for LINE, SINE, RETROPOSON
void filterLINE(Cluster *clt)
{
    // For all cluster
    if ((hasPolyA(clt->flag) || isATRich(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) && hasTSD(clt->flag))
        clt->flag |= CLT_PASS;

    // For cluster locates in normal region
    if (!isSelfToSelf(clt->flag)) {
        if ((hasPolyA(clt->flag) || isATRich(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
            clt->flag |= CLT_PASS;

        if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) && hasSingleTE(clt->flag) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
            clt->flag |= CLT_PASS;
    }
}

// Perform post-filtering for LTR
void filterLTR(Cluster *clt)
{
    // For all cluster
    if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) && hasTSD(clt->flag))
        clt->flag |= CLT_PASS;

    if (isSoloLtr(clt->flag) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) && hasTSD(clt->flag))
        clt->flag |= CLT_PASS;

    // For cluster locates in normal region
    if (!isSelfToSelf(clt->flag)) {
        if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
            clt->flag |= CLT_PASS;

        if (isSoloLtr(clt->flag) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
            clt->flag |= CLT_PASS;
    }
}

// Perform post-filtering for DNA
void filterDNA(Cluster *clt)
{
    // For all cluster
    if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) && hasTSD(clt->flag))
        clt->flag |= CLT_PASS;

    // For cluster locates in normal region
    if (!isSelfToSelf(clt->flag)) {
        if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
            clt->flag |= CLT_PASS;
    }
}