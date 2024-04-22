
#include "baseline/basics.h"

void *myalloc (size_t n)

   { void *p;
     if (n == 0) return NULL;
     p = malloc(n);
     if (p == NULL) 
        { fprintf(stderr,"Error: malloc of %li bytes returned null\n",n);
	  exit(1);
	}
     return p;
   }

void *mycalloc (size_t n, size_t b)

{ void *p;
    if (n == 0) return NULL;
    p = calloc(n, b);
    if (p == NULL)
    { fprintf(stderr,"Error: calloc of %li bytes returned null\n",n*b);
        exit(1);
    }
    return p;
}

void *myrealloc (void *p, size_t n)

   { if (p == NULL) return myalloc(n);
     if (n == 0) return NULL;
     p = realloc(p,n);
     if (p == NULL)
        { fprintf(stderr,"Error: realloc of %li bytes returned null\n",n);
	  exit(1);
	}
     return p;
   }


void myfree (void *p)

  { if (p != NULL) free(p);
  }

uint numbits (uint64_t n)

   { uint bits = 0;
     while (n)
        { n >>= 1; bits++; }
     return bits ? bits : 1;
   }

        // accesses A[i] from b-bit packed array A

extern inline uint64_t access64 (uint64_t *A, uint64_t pos, uint b)

   { uint64_t bstart = pos*b;
     uint64_t pstart = bstart / w;
     uint64_t offs = bstart % w;
     uint64_t val = A[pstart] >> offs;
     if (offs+b > w) val |= A[pstart+1] << (w-offs);
     return val & ((((uint64_t)1) << b)-1);
   }


        // writes A[i]=val in b-bit packed array A

extern inline void write64 (uint64_t *A, uint64_t pos, uint64_t val, uint b)

   { uint64_t bstart = pos*b;
     uint64_t pstart = bstart / w;
     uint64_t offs = bstart % w;
     uint64_t cell;
     cell = A[pstart] & ~(((((uint64_t)1) << b)-1) << offs);
     A[pstart] = cell | (val << offs);
     if (offs+b > w)
	{ cell = A[pstart+1] & ~((((uint64_t)1) << (offs+b-w))-1);
          A[pstart+1] = cell | (val >> (w-offs));
	}
   }

