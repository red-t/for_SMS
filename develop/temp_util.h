#ifndef TEMP_UTIL_H
#define TEMP_UTIL_H

#include "htslib/sam.h"


/*! @function
 @abstract  Get whether the query is secondary alignment
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
#define bam_is_second(b) (((b)->core.flag&BAM_FSECONDARY) != 0)


#endif // TEMP_UTIL_H