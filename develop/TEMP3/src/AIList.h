//=====================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  09/5/2018
//modified by Zhongren Hu	14/3/2023
//-------------------------------------------------------------------------------------
#ifndef __AILIST_H__
#define __AILIST_H__
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
//-------------------------------------------------------------------------------------
#include "khash.h"
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAXC 10							//max number of components
//-------------------------------------------------------------------------------------
typedef struct {
    int32_t start;      				// region start: 0-based
    int32_t end;    					// region end: not inclusive
                                        // TO DO: should add a new feature value, representing the TE family
} gdata_t;

typedef struct{
	int64_t nr, mr;						// number of regions; max number
	gdata_t *glist;						// list of regions
	int nc, lenC[MAXC], idxC[MAXC];		// number of components; lengths of components; start indices of components 
	int32_t *maxE;						// list of maximum end (maxE)
} ctg_t;

typedef struct {
	ctg_t *ctg;        					// list of contigs (only 1 contig)
} ailist_t;

//-------------------------------------------------------------------------------------

/// Initialize ailist_t
ailist_t *ailist_init(void);


/// Free ailist data.
void ailist_destroy(ailist_t *ail);


/// Add a interval into AIList.
/*!
 * @param s   start of the interval
 * @param e   end of the interval
 */
void ailist_add(ailist_t *ail, int32_t s, int32_t e);


/// Add intervals from BED-like file.
/*!
 * @param fn   filename of the BED-like file
 * @param chr  name of the specified chromosome, only intervals
 *             on this chromosome will be added
 */
void readBED(ailist_t *ail, const char* fn, const char* chr);


/// Construct ailist: decomposition and augmentation.
/*!
 * @param cLen minimum coverage length, default=20
 */
void ailist_construct(ailist_t *ail, int cLen);


/// Binary search
/*!
 * @param As    pointer to the list of intervals
 * @param idxS  start index of the search space
 * @param idxE  end index of the search space
 * @param qe    query end for searching index of the search space
 * @return      index of the first item satisfying .start<qe from right
 */
int32_t bSearch(gdata_t* As, int32_t idxS, int32_t idxE, int32_t qe);


/// Compute the minimum distance between query position to the nearest interval
/*!
 * @param rpos  reference position used as query
 * @param flank flank size used for extending query
 * @param n     number of intervals overlapped with [rpos-flank, rpos+flank)
 * @param d     minimum distance, will be computed when calling this function
 */
void query_dist_p(ailist_t *ail, int32_t rpos, int32_t flank, int32_t *n, int32_t *d);


/// Compute the minimum distance between query interval to the nearest interval
/*!
 * @param st    start position of query interval
 * @param ed    end position of query interval
 * @param flank flank size used for extending query interval
 * @param d     minimum distance, will be computed when calling this function
 * @return      number of intervals overlapped with [st-flank, ed+flank)
 */
void query_dist_c(ailist_t *ail, int32_t st, int32_t ed, int32_t flank, int32_t *n, int32_t *d);

/*********************
 * Convenient macros *
 *********************/
#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	}while (0) 

#endif