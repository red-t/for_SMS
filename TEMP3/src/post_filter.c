#include "post_filter.h"

// Perform post-filtering for different TE class
void postFilter(Cluster *clt)
{
    uint32_t CLT_CLASS_MASK = 0xf8000;
    switch (clt->flag & CLT_CLASS_MASK)
    {
    case CLT_LINE:
    case CLT_SINE:
    case CLT_RETROPOSON:
        filterLINE(clt);
        break;
    case CLT_LTR:
        filterLTR(clt);
        break;
    case CLT_DNA:
        filterDNA(clt);
        break;
    default:
        filterDNA(clt);
        break;
    }
}

// Perform post-filtering for LINE, SINE, RETROPOSON
void filterLINE(Cluster *clt)
{
    // For all cluster
    if ((hasPolyA(clt->flag) || isATRich(clt->flag)) && (isLeftNearEnd(clt->flag) && isRightNearEnd(clt->flag)))
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