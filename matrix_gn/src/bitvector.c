
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bitvector.h"

extern inline uint popcount (uint64_t y)

    { y -= ((y >> 1) & 0x5555555555555555ull);
      y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
      return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
    }

	// creates an empty bitvector of n bits

bitvector bitsCreate (uint64_t n)

    { bitvector B;
      B = (bitvector)malloc(sizeof(struct s_bitvector));
      B->size = n;
      if (n == 0) B->data = NULL;
      else B->data = (uint64_t*)malloc(((n+w-1)/w)*sizeof(uint64_t));
      B->owned = 1;
      B->k = 0;
      B->S = NULL;
      return B;
    }

	// converts a bit array into a bitvector of n bits
        // data is pointed to and will be freed if own


bitvector bitsCreateFrom (void *data, uint64_t n, uint own)

    { bitvector B;
      B = (bitvector)malloc(sizeof(struct s_bitvector));
      B->size = n;
      if (n == 0) B->data = NULL;
      else B->data = (uint64_t*)data;
      B->owned = own;
      B->k = 0;
      B->S = NULL;
      return B;
    }

	// destroys B

void bitsDestroy (bitvector B)

    { if (B != NULL) 
         { if (B->owned && (B->data != NULL)) free(B->data);
           if (B->S != NULL) free(B->S);
	   B->data = NULL; B->S = NULL; // for safety
      	   free(B);
	 }
    }

        // creates a new copy of B, with its own data and same preprocessing

bitvector bitsCopy (bitvector B)

   { uint64_t siz;
     bitvector C = (bitvector)malloc(sizeof(struct s_bitvector));
     *C = *B;
     if (B->data != NULL)
        { siz = ((C->size+w-1)/w)*sizeof(uint64_t);
	  C->data = (uint64_t*)malloc(siz);
          memcpy(C->data,B->data,siz);
	}
     if (B->S != NULL)
        { siz = ((B->size+C->k*w-1)/(C->k*w))*sizeof(uint64_t);
          C->S = (uint64_t*)malloc(siz);
          memcpy(C->S,B->S,siz);
	}
     C->owned = 1;
     return C;
   }

        // writes B to file, which must be opened for writing 

void bitsSave (bitvector B, FILE *file)

   { fwrite (&B->size,sizeof(uint64_t),1,file);
     if (B->size != 0)
        fwrite (B->data,sizeof(uint64_t),(B->size+w-1)/w,file);
   }

        // loads bitvector from file, which must be opened for reading
	// must be later preprocessed for rank

bitvector bitsLoad (FILE *file)

    { bitvector B;
      B = (bitvector)malloc(sizeof(struct s_bitvector));
      fread (&B->size,sizeof(uint64_t),1,file);
      if (B->size == 0) B->data = NULL;
      else { B->data = (uint64_t*)malloc(((B->size+w-1)/w)*sizeof(uint64_t));
	     fread (B->data,sizeof(uint64_t),(B->size+w-1)/w,file);
	   }
      B->owned = 1;
      B->k = 0;
      B->S = NULL;
      return B;
    }

	// bitvector size in w-bit words

uint64_t bitsSpace (bitvector B)

    { uint64_t space = sizeof(struct s_bitvector)*8/w;
      if (B == NULL) return 0;
      if (B->data != NULL) space += (B->size+w-1)/w;
      if (B->S != NULL) space += (B->size+B->k*w-1)/(B->k*w);
      return space;
    }

        // gives bit data

extern inline uint64_t *bitsData (bitvector B)

    { return B->data;
    }

        // gives bit length

extern inline uint64_t bitsLength (bitvector B)

    { return B->size;
    }

	// sets value for B[i]= (v != 0), assumes i is right

extern inline void bitsWriteA (uint64_t *B, uint64_t i, uint v)

    { if (v == 0) B[i/w] &= ~(((uint64_t)1) << (i%w));
      else B[i/w] |= ((uint64_t)1) << (i%w);
    }

extern inline void bitsWrite (bitvector B, uint64_t i, uint v)

    { bitsWriteA (B->data,i,v);
    }

	// access B[i], assumes i is right

extern inline uint bitsAccessA (uint64_t *B, uint64_t i)

    { return (B[i/w] >> (i%w)) & 1;
    }

extern inline uint bitsAccess (bitvector B, uint64_t i)

    { return bitsAccessA(B->data,i);
    }

	// preprocesses for rank, with parameter k, no writes after this

void bitsRankPreprocess (bitvector B, uint64_t k)

    { uint64_t i,n;
      uint64_t acc;
      n = B->size;
      if (n == 0) return;
      if (B->S != NULL) free(B->S); // a previous preprocessing
      B->S = (uint64_t*)malloc(((n+k*w-1)/(k*w))*sizeof(uint64_t));
      B->k = k;
      acc = 0;
      for (i=0;i<(n+w-1)/w;i++)
	  { if (i%k == 0) B->S[i/k] = acc;
	    acc += popcount(B->data[i]);
	  }
    } 

	// computes rank(B,i), zero-based, assumes i is right

extern inline uint64_t bitsRank (bitvector B, uint64_t i)

    { uint64_t b,sb;
      uint64_t rank;
      sb = i/(B->k*w);
      rank = B->S[sb];
      sb *= B->k;
      for (b=sb;b<i/w;b++) rank += popcount(B->data[b]);
      return rank + popcount(B->data[b] & (((uint64_t)~0) >> (w-1-(i%w))));
    }

