
#ifndef INCLUDEDbasics
#define INCLUDEDbasics

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char byte;

typedef uint32_t uint;

#define w (8*sizeof(uint64_t))
#define mmin(x,y) (((x)<(y))?(x):(y))
#define mmax(x,y) (((x)>(y))?(x):(y))

void *myalloc (size_t n);
void *myrealloc (void *p, size_t n);
void myfree (void *p);

	// number of bits needed to represent n, gives 1 for n=0
uint numbits (uint n);

#endif
