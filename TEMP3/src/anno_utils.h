#include <stdint.h>

/******************
 *** Data Types ***
 ******************/

typedef struct Anno
{
    int     idx;
    int     refStart;
    int     refEnd;
    int     tid;
    uint8_t strand;
    int     queryStart;
    int     queryEnd;
} __attribute__((packed)) Anno;