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
    int         idx;
    int         cltTid;
    int         cltIdx;
    int         queryStart;
    int         queryEnd;
    uint8_t     strand;
    int         tid;
    int         refStart;
    int         refEnd;
    uint32_t    flag;
} __attribute__((packed)) Anno;

/// @brief Data container for polyA
typedef struct PolyA
{
    int idx;
    int cltTid;
    int cltIdx;
    int leftAnnoStart;
    int leftIdx;
    int rightAnnoEnd;
    int rightIdx;
    int isA;
    int seqLen;
} PolyA;


/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/
#define isRevAnno(anno) ((anno).strand == 1)

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fillAnnoArray(Cluster *clt, Anno *annoArray, int idx);

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *clt, Anno *annoArray, int numAnno, PolyA *polyA);

/// @brief Find and record single polyA/polyT
int setPolyA(char *flankSeq, Anno *annoArray, Cluster *clt, int numAnno, PolyA *polyA);

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *clt, PolyA *polyA, Anno *annoArray, int numAnno);

/// @brief Adjust annotation position
void adjustAnno(Anno *annoArray, int numAnno, int leftDelta);

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annoTsd(Cluster *clt, Anno *annoArray, int numAnno);

/// @brief Find TSD and refine breakpoint
int setTsd(Cluster *clt, int localStart, int leftEnd, int rightStart);

/// @brief Set ins-seq structure based on annotations
void setInsStruc(Cluster *clt, Anno *annoArray, int numAnno, uint32_t *classArray, int *sizeArray, int *ltrArray);

/// @brief Compare function for sorting annotations
int compare(const void *a, const void *b);

/// @brief Check whether the ins-seq contains large gap
void checkGap(Anno *annoArray, int numAnno, Cluster *clt, int *sizeArray);

/// @brief Check whether the ins-seq contains valid polyA
void checkPolyA(Anno *annoArray, int numAnno, Cluster *clt);

/// @brief Check whether the ins-seq contains complete ends
void checkEnd(Anno *annoArray, int numAnno, Cluster *clt);

/// @brief Check which TE class the insertion belongs to
void checkTEClass(Anno *annoArray, int numAnno, Cluster *clt, uint32_t *classArray);

/// @brief Check whether left-/right- assm-flank-seq contains valid polyT/A
void checkFlankPolyA(Anno *annoArray, int numAnno, Cluster *clt);

/// @brief Search polyT/polyA in left-/right- assm-flank-seq sequence
int searchFlankPolyA(char *flankSeq, int isA, int seqLen);

/// @brief Check whether the insertion is SOLO LTR
void checkSoloLtr(Anno *annoArr, int numAnno, Cluster *clt, int *sizeArr, int *ltrArr);


/**********************
 *** Annotation I/O ***
 **********************/

/// @brief Output formated annotation records
void outputAnno(Anno *annoArray, int numAnno, int startIdx, const char *teFn);

/// @brief Change single annotation record into specified format
void formatSingleAnno(Anno anno, char *queryTmp, char *refTmp, faidx_t *teFa, int *teTable, int *strandFlag);

/// @brief Write annotation for a single cluster
void writeSingleCltAnno(int strandFlag, int numTe, int *teTable, faidx_t *teFa, char *queryStr, char *refStr, FILE *fp, Anno anno);

/// @brief Get cluster strand
char getCltStrand(int strandFlag);

/// @brief Generate a string which represents cluster TE-class
char *getCltClass(int numTe, int *teTable, faidx_t *teFa);

#endif // ANNO_UTILS_H