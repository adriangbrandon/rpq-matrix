
#ifndef INCLUDEDbitvector
#define INCLUDEDbitvector

	// supports bitvectors of size up to 2^64-1

#include <stdint.h>

#include <stdio.h>

typedef unsigned char byte;
typedef unsigned uint;
#define w (8*sizeof(uint64_t))
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))

typedef struct s_bitvector {
    uint64_t size; // number of bits
    uint64_t* data; // the bits
    int owned; // data is owned (to be freed)
    uint64_t k; // superblock length is 64 k
    uint64_t *S; // superblocks
    } *bitvector;

	// creates an empty bitvector of n bits
bitvector bitsCreate (uint64_t n);

	// converts a bit array into a bitvector of n bits
	// data is pointed to and will be freed if own
bitvector bitsCreateFrom (void *data, uint64_t n, uint own);

	// destroys B, frees data if it was owned
void bitsDestroy (bitvector B);

        // creates a new copy of B, with its own data and same preprocessing
bitvector bitsCopy (bitvector B);

	// writes B to file, which must be opened for writing
void bitsSave (bitvector B, FILE *file);

	// loads bitvector from file, which must be opened for reading
        // must be later preprocessed for rank
bitvector bitsLoad (FILE *file);

	// gives space of bitvector in w-bit words
uint64_t bitsSpace (bitvector B);

	// gives bit data
extern inline uint64_t *bitsData (bitvector B);

	// gives bit length
extern inline uint64_t bitsLength (bitvector B);

	// sets value for B[i]= (v != 0), assumes i is right
extern inline void bitsWrite (bitvector B, uint64_t i, uint v);

	// a more basic version on arrays 
extern inline void bitsWriteA (uint64_t *B, uint64_t i, uint v);

	// access B[i], assumes i is right
extern inline uint bitsAccess (bitvector B, uint64_t i);

	// a more basic version on arrays 
extern inline uint bitsAccessA (uint64_t *B, uint64_t i);

        // counts # of 1s in y

extern inline uint popcount (uint64_t y);

	// preprocesses for rank, with parameter k, no writes after this
	// frees previous S if it exists
void bitsRankPreprocess (bitvector B, uint64_t k);

	// computes rank(B,i), zero-based, assumes i is right
extern inline uint64_t bitsRank (bitvector B, uint64_t i);

#endif
