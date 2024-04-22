
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
void *mycalloc (size_t n, size_t b);
void *myrealloc (void *p, size_t n);
void myfree (void *p);

	// number of bits needed to represent n, gives 1 for n=0
uint numbits (uint64_t n);

	// w-words needed to allocate for n b-bit values
#define packedwords(n,b) (((n)*((uint64_t)(b))+w-1)/w)

	// accesses A[i] from b-bit packed array A
extern inline uint64_t access64 (uint64_t *A, uint64_t pos, uint b);

	// writes A[i]=val in b-bit packed array A
extern inline void write64 (uint64_t *A, uint64_t pos, uint64_t val, uint b);

#endif
