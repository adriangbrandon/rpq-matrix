
#include "basics.h"

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

uint numbits (uint n)

   { uint bits = 0;
     while (n)
        { n = n>>1; bits++; }
     return bits ? bits : 1;
   }

