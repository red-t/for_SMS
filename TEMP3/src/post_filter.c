#include "post_filter.h"

// Perform post-filtering for different TE class
void postFilter(Cluster *clt)
{
    if ((clt->flag & (CLT_LINE|CLT_SINE|CLT_RETROPOSON)) != 0)
        filterLINE(clt);
    
    if ((clt->flag & (CLT_LTR|CLT_DNA)) != 0)
        filterLTR(clt);
}

// Perform post-filtering for LINE, SINE, RETROPOSON
void filterLINE(Cluster *clt)
{
    // For all cluster
    if (hasPolyA(clt->flag) || isATRich(clt->flag))
        if (hasSingleTE(clt->flag) || (hasFull5P(clt->flag) && hasFull3P(clt->flag)) || (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) || (hasTSD(clt->flag)))
            clt->flag |= CLT_PASS;
    
    if (hasSingleTE(clt->flag))
        if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) || (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) || (hasTSD(clt->flag)))
            clt->flag |= CLT_PASS;

    if (hasFull5P(clt->flag) && hasFull3P(clt->flag))
        if ((isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) || (hasTSD(clt->flag)))
            clt->flag |= CLT_PASS;

    if (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag))
        if (hasTSD(clt->flag))
            clt->flag |= CLT_PASS;

    // For cluster locates in normal region
    if (!isSelfToSelf(clt->flag)) {
        if (hasPolyA(clt->flag) || isATRich(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasSingleTE(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasFull5P(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;

        if (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;
        
        if (hasTSD(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasSingleTE(clt->flag) && hasFull5P(clt->flag) && hasUnknown3P(clt->flag))
            clt->flag |= CLT_PASS;
    }
}

// Perform post-filtering for LTR, DNA
void filterLTR(Cluster *clt)
{
    // For all cluster
    if (hasSingleTE(clt->flag))
        if ((hasFull5P(clt->flag) && hasFull3P(clt->flag)) || (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) || (hasTSD(clt->flag)))
            clt->flag |= CLT_PASS;

    if (hasFull5P(clt->flag) && hasFull3P(clt->flag))
        if ((isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)) || (hasTSD(clt->flag)))
            clt->flag |= CLT_PASS;

    if (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag))
        if (hasTSD(clt->flag))
            clt->flag |= CLT_PASS;

    // For cluster locates in normal region
    if (!isSelfToSelf(clt->flag))
    {
        if (hasFull5P(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;

        if (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasTSD(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasSingleTE(clt->flag) && hasFull5P(clt->flag))
            clt->flag |= CLT_PASS;

        if (hasSingleTE(clt->flag) && hasFull3P(clt->flag))
            clt->flag |= CLT_PASS;
    }
}