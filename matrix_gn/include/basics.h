
#ifndef INCLUDEDbasics
#define INCLUDEDbasics

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char byte;
typedef unsigned uint;

#define w (8*sizeof(uint64_t))
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))

void *myalloc (size_t n);
void *myrealloc (void *p, size_t n);
void myfree (void *p);

#endif
