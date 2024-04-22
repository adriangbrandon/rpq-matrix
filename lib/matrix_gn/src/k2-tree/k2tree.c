
#include <stdlib.h>
#include <string.h>

#include "k2-tree/k2tree.h"

#define rankK 4  // parameter for bit rank

uint mapId[] = { 0, 1, 2, 3 };
uint mapTr[] = { 0, 2, 1, 3 };
static uint remapTr[] =
		{ 0, 1, 4, 5,  2, 3, 6, 7,  8, 9, 12, 13,  10, 11, 14, 15 };

	// creates k2tree from n coords (2n ints x,y) that need nbits bits
	// reorders coords array

typedef struct {
   uint64_t *coords1, *coords2;
   uint64_t n;
   uint level;
   } queue;

typedef struct {
   uint32_t *coords1, *coords2;
   uint64_t n;
   uint level;
   } queue32;

static queue *Q = NULL;
static queue32 *Q32 = NULL;
static uint64_t head,tail,size;

static uint64_t *B = NULL;

	// computes levels of the k2tree

static void k2computeLevels (k2tree T)

   { int l;
     uint64_t p;
     T->levels = (uint64_t*)myalloc(T->nlevels * sizeof(uint64_t));
     p = 1;
     for (l=T->nlevels-1;l>=0;l--)
	 { T->levels[l] = p;
	   p = bitsRank(T->B,4*p-1)+1;
	 }
     T->elems = p-T->levels[0];
   }

static void k2write (uint level, uint64_t ptr, uint64_t n,
		     uint64_t *coords1, uint64_t *coords2)
   { uint64_t cnt[5];
     uint64_t i,prev;
     uint v;
     for (v=0;v<=4;v++) cnt[v] = 0;
     for (i=0;i<2*n;i+=2)
	 { v = 2 * ((coords1[i] >> level) & 1) 
	       + ((coords1[i+1] >> level) & 1);
	   cnt[v+1]++;
	 }
     for (v=1;v<=4;v++)
         bitsWriteA(B,ptr++,cnt[v]);
     if (level)
        { cnt[2] += cnt[1]; cnt[3] += cnt[2];
          for (i=0;i<2*n;i+=2)
	      { v = 2 * ((coords1[i] >> level) & 1) 
		    + ((coords1[i+1] >> level) & 1);
	        coords2[2*cnt[v]] = coords1[i];
	        coords2[2*cnt[v]+1] = coords1[i+1];
	        cnt[v]++;
	      }
	  prev = 0;
          for (v=0;v<4;v++)
	      { n = cnt[v]-prev;
	        if (n) 
	           { Q[tail].level = level-1;
	             Q[tail].n = n;
	             Q[tail].coords1 = coords2+2*prev;
	             Q[tail].coords2 = coords1+2*prev;
	             tail = (tail+1)%size;
		     prev = cnt[v];
	           }
	      }
	}
   }

static void k2write32 (uint level, uint64_t ptr, uint64_t n,
		       uint32_t *coords1, uint32_t *coords2)
   { uint64_t cnt[5];
     uint64_t i,prev;
     uint v;
     for (v=0;v<=4;v++) cnt[v] = 0;
     for (i=0;i<2*n;i+=2)
	 { v = 2 * ((coords1[i] >> level) & 1) 
	       + ((coords1[i+1] >> level) & 1);
	   cnt[v+1]++;
	 }
     for (v=1;v<=4;v++)
         bitsWriteA(B,ptr++,cnt[v]);
     if (level)
        { cnt[2] += cnt[1]; cnt[3] += cnt[2];
          for (i=0;i<2*n;i+=2)
	      { v = 2 * ((coords1[i] >> level) & 1) 
		    + ((coords1[i+1] >> level) & 1);
	        coords2[2*cnt[v]] = coords1[i];
	        coords2[2*cnt[v]+1] = coords1[i+1];
	        cnt[v]++;
	      }
	  prev = 0;
          for (v=0;v<4;v++)
	      { n = cnt[v]-prev;
	        if (n) 
	           { Q32[tail].level = level-1;
	             Q32[tail].n = n;
	             Q32[tail].coords1 = coords2+2*prev;
	             Q32[tail].coords2 = coords1+2*prev;
	             tail = (tail+1)%size;
		     prev = cnt[v];
	           }
	      }
	}
   }

k2tree k2create (uint64_t n, uint nbits, uint64_t *coords)

   { k2tree T;
     uint64_t i,len;
     uint64_t *coords2;

     T = (k2tree)myalloc(sizeof(struct s_k2tree));
     T->nlevels = nbits;
     T->levels = (uint64_t*)myalloc(nbits * sizeof(uint64_t));
     B = (uint64_t*)myalloc(((4*(1+n*nbits)+w-1)/w)*sizeof(uint64_t));
     size = n+1;
     Q = (queue*)myalloc(size*sizeof(queue));
     coords2 = (uint64_t*)myalloc(2*n*sizeof(uint64_t));
     Q[0].coords1 = coords; Q[0].coords2 = coords2;
     Q[0].n = n; Q[0].level = nbits-1;
     head = 0; tail = 1; 
     len = 0;
     while (head != tail)
	{ k2write (Q[head].level,len,Q[head].n,
		   Q[head].coords1,Q[head].coords2);
	  head = (head+1)%size;
	  len += 4;
	}
     myfree (coords2);
     myfree(Q); Q = NULL; // safety
     B = (uint64_t*)myrealloc(B,((len+w-1)/w)*sizeof(uint64_t));
     T->B = bitsCreateFrom(B,len,1);
     B = NULL; // safety
     bitsRankPreprocess(T->B,rankK);
     k2computeLevels (T);
     return T;
   }

k2tree k2create32 (uint64_t n, uint nbits, uint32_t *coords)

   { k2tree T;
     uint64_t i,len;
     uint32_t *coords2;

     T = (k2tree)myalloc(sizeof(struct s_k2tree));
     T->nlevels = nbits;
     T->levels = (uint64_t*)myalloc(nbits * sizeof(uint64_t));
     B = (uint64_t*)myalloc(((4*(1+n*nbits)+w-1)/w)*sizeof(uint64_t));
     size = n+1;
     Q32 = (queue32*)myalloc(size*sizeof(queue32));
     coords2 = (uint32_t*)myalloc(2*n*sizeof(uint32_t));
     Q32[0].coords1 = coords; Q32[0].coords2 = coords2;
     Q32[0].n = n; Q32[0].level = nbits-1;
     head = 0; tail = 1; 
     len = 0;
     while (head != tail)
	{ k2write32 (Q32[head].level,len,Q32[head].n,
		     Q32[head].coords1,Q32[head].coords2);
	  head = (head+1)%size;
	  len += 4;
	}
     myfree (coords2);
     myfree(Q32); Q32 = NULL; // safety
     B = (uint64_t*)myrealloc(B,((len+w-1)/w)*sizeof(uint64_t));
     T->B = bitsCreateFrom(B,len,1);
     B = NULL; // safety
     bitsRankPreprocess(T->B,rankK);
     k2computeLevels (T);
     return T;
   }

	// creates k2tree from bits, of length len, nlevels levels
        // if own, will free bitvector when destroyed

k2tree k2createFrom (uint nlevels, uint64_t len, void *bits, uint own)

   { k2tree T = (k2tree)myalloc(sizeof(struct s_k2tree));
     T->nlevels = nlevels;
     T->B = bitsCreateFrom(bits,len,own);
     bitsRankPreprocess(T->B,rankK);
     k2computeLevels (T);
     return T;
   }

	// destroys k2tree

void k2destroy (k2tree T)

   { if (T == NULL) return;
     myfree(T->levels);
     bitsDestroy (T->B);
     myfree(T);
   }

        // creates a new copy of T, with its own data

k2tree k2copy (k2tree T)

   { k2tree C = (k2tree)myalloc(sizeof(struct s_k2tree));
     *C = *T;
     C->levels = (uint64_t*)myalloc(C->nlevels*sizeof(uint64_t));
       memcpy(C->levels,T->levels,C->nlevels*sizeof(uint64_t));
     C->B = bitsCopy(T->B);
     return C;
   }

        // writes T to file, which must be opened for writing

void k2save (k2tree T, FILE *file)

   { fwrite (&T->nlevels,sizeof(uint),1,file);
     bitsSave (T->B,file);
   }

        // loads k2tree from file, which must be opened for reading

k2tree k2load (FILE *file)

    { k2tree T;
      T = (k2tree)myalloc(sizeof(struct s_k2tree));
      fread (&T->nlevels,sizeof(uint),1,file);
      T->B = bitsLoad(file);
      bitsRankPreprocess(T->B,rankK);
      k2computeLevels(T);
      return T;
    }

	// space of the k2tree, in w-bit words

uint64_t k2space (k2tree T)

   { uint64_t space = sizeof(struct s_k2tree)*8/w;
     if (T == NULL) return 0;
     if (T->levels != NULL) space += T->nlevels;
     if (T->B != NULL) space += bitsSpace(T->B);
     return space;
   }

	// levels of the k2tree

extern inline uint k2levels (k2tree T)

   { return T->nlevels;
   }

	// elements in the k2tree

extern inline uint64_t k2elems (k2tree T)

   { return T->elems;
   }

	// bitvector of the k2tree

extern inline bitvector k2bits (k2tree T)

   { return T->B;
   }

	// root of the k2tree

extern inline k2node k2root (k2tree T)

   { return 0;
   }

	// u has i-th child, assumes u is right and i = 0..3

extern inline uint k2hasChild (k2tree T, k2node u, uint i)

   { return bitsAccess(T->B,4*u+i);
   }

	// returns the 4-bit signature of node u

extern inline uint k2sigNode (k2tree T, k2node u)

   { uint64_t i = 4*u;
     return (T->B->data[i/w] >> (i%w)) & 0xF;
   }

	// i-th child of u, assumes u is right, i = 0..3, and child exists

extern inline k2node k2child (k2tree T, k2node u, uint i)

   { return bitsRank(T->B,4*u+i);
   }

        // fills a child vector for the 4 children, returns signature
	// assumes 4 divides w

extern inline uint k2fillChildren (k2tree T, k2node u, k2node *child)

   { uint64_t i = 4*u;
     uint64_t r = i ? bitsRank(T->B,i-1) : 0;
     uint s = (T->B->data[i/w] >> (i%w)) & 0xF;
     uint j,v,sig = s;
     for (j=0;j<4;j++)
	 { v = s & 1;
	   r += v;
	   s >>= 1;
	   child[j] = r;
	 }
     return sig;
   }

	// version that allows mapped children

extern inline uint k2fillMappedChildren (k2tree T, k2node u, uint *map,
			                 k2node *child)

   { uint64_t i = 4*u;
     uint64_t r = i ? bitsRank(T->B,i-1) : 0;
     uint s = (T->B->data[i/w] >> (i%w)) & 0xF;
     uint j,v,sig;
     sig = (map == mapId) ? s : remapTr[s];
     for (j=0;j<4;j++)
	 { v = s & 1;
	   r += v;
	   s >>= 1;
		// I'm applying the inverse perm here, but for transp is ok
	   child[map[j]] = r;
	 }
     return sig;
   }

        // recovers all the cells in [r1..r2] x [c1..c2]
        // writes 2n integers in buffer, which must be preallocated
        // to 2*elems of uint64_t (worst case)
        // returns number of elements found. just counts if buffer is NULL

static uint64_t collect (k2tree T, k2node u, uint level, 
		         uint64_t r1, uint64_t r2, uint64_t c1, uint64_t c2, 
		         uint64_t roffs, uint64_t coffs, 
			 uint64_t *buffer, uint64_t count, uint cr)

   { uint64_t lim;
     if (level==0)
	{ if (buffer != NULL)
	     { buffer[2*count+cr] = roffs; 
	       buffer[2*count+1-cr] = coffs; 
	     }
	  return count+1;
	}
     level--;
     lim = ((uint64_t)1)<<level;
     if ((r1 < lim) && (c1 < lim) && k2hasChild(T,u,0))
	{ count = collect(T,k2child(T,u,0),level,
			   r1,mmin(r2,lim-1),c1,mmin(c2,lim-1),
			   roffs,coffs,buffer,count,cr);
	}
     if ((r1 < lim) && (c2 >= lim) && k2hasChild(T,u,1))
	{ count = collect(T,k2child(T,u,1),level,
			   r1,mmin(r2,lim-1), mmax(c1,lim)-lim,c2-lim,
			   roffs,coffs+lim,buffer,count,cr);
	}
     if ((r2 >= lim) && (c1 < lim) && k2hasChild(T,u,2))
	{ count = collect(T,k2child(T,u,2),level,mmax(r1,lim)-lim,r2-lim,c1,mmin(c2,lim-1),
			   roffs+lim,coffs,buffer,count,cr);
	}
     if ((r2 >= lim) && (c2 >= lim) && k2hasChild(T,u,3))
	{ count = collect(T,k2child(T,u,3),level,
			   mmax(r1,lim)-lim,r2-lim,mmax(c1,lim)-lim,c2-lim,
			   roffs+lim,coffs+lim,buffer,count,cr);
	}
     return count;
   }

uint64_t k2collect (k2tree T, uint64_t r1, uint64_t r2, 
		    uint64_t c1, uint64_t c2, uint64_t *buffer, uint cr)

   { return collect(T,k2root(T),k2levels(T),r1,r2,c1,c2,0,0,buffer,0,cr);
   }

static uint64_t collect32 (k2tree T, k2node u, uint level, 
		         uint32_t r1, uint32_t r2, uint32_t c1, uint32_t c2, 
		         uint32_t roffs, uint32_t coffs, 
			 uint32_t *buffer, uint64_t count, uint cr)

   { uint64_t lim;
     if (level==0)
	{ if (buffer != NULL)
	     { buffer[2*count+cr] = roffs; 
	       buffer[2*count+1-cr] = coffs; 
	     }
	  return count+1;
	}
     level--;
     lim = ((uint64_t)1)<<level;
     if ((r1 < lim) && (c1 < lim) && k2hasChild(T,u,0))
	{ count = collect32(T,k2child(T,u,0),level,
			   r1,mmin(r2,lim-1),c1,mmin(c2,lim-1),
			   roffs,coffs,buffer,count,cr);
	}
     if ((r1 < lim) && (c2 >= lim) && k2hasChild(T,u,1))
	{ count = collect32(T,k2child(T,u,1),level,
			   r1,mmin(r2,lim-1),mmax(c1,lim)-lim,c2-lim,
			   roffs,coffs+lim,buffer,count,cr);
	}
     if ((r2 >= lim) && (c1 < lim) && k2hasChild(T,u,2))
	{ count = collect32(T,k2child(T,u,2),level,
			   mmax(r1,lim)-lim,r2-lim,c1,mmin(c2,lim-1),
			   roffs+lim,coffs,buffer,count,cr);
	}
     if ((r2 >= lim) && (c2 >= lim) && k2hasChild(T,u,3))
	{ count = collect32(T,k2child(T,u,3),level,
			   mmax(r1,lim)-lim,r2-lim,mmax(c1,lim)-lim,c2-lim,
			   roffs+lim,coffs+lim,buffer,count,cr);
	}
     return count;
   }

uint64_t k2collect32 (k2tree T, uint32_t r1, uint32_t r2, 
		    uint32_t c1, uint32_t c2, uint32_t *buffer, uint cr)

   { return collect32(T,k2root(T),k2levels(T),r1,r2,c1,c2,0,0,buffer,0,cr);
   }

	// merges the bit arrays of two k2trees, writes its bit length
	// in *len and number of elements in *telems. nodes per level are
	// written in array *levels if not null, levels[0] is #leaves

typedef unsigned short queue2; // just one short as this is space-critical
			// limits levels to 64, which should be more than ok

#define addnumc(x,o) x += (o)<<8
#define getnumc(x) (1+((x)>>8))
#define gethas(x) (((x)>>6)&0x03)
#define getdist(x) ((x)&0x3F)
#define hasdist(h,n,d) ((((n)-1) << 8) | ((h) << 6) | (d))

// copies v from tree[ptr..] in bits

static inline void addWord (uint64_t *tree, uint64_t ptr, uint64_t v)

   { uint m = ptr % w;
     if (m)
        { tree[ptr/w] = (tree[ptr/w] & ((((uint64_t)1) << m) - 1)) | (v << m);
          tree[ptr/w+1] = (tree[ptr/w+1] & ((~(uint64_t)0) << m)) | (v>>(w-m));
	}
     else tree[ptr/w] = v;
   }

uint64_t *k2merge (uint64_t *treeA, uint64_t lenA, 
		   uint64_t* treeB, uint64_t lenB, uint level, 
		   uint64_t *len, uint64_t *telems, uint64_t *levels)

   { uint64_t *treeM;
     uint64_t ptr,ptrA,ptrB;
     uint v,sig,sigA,sigB;
     uint64_t sigs;
     uint64_t elems,size;
     queue2 *Q;
     uint64_t head,tail,prev;
     uint o,n;

     treeM = (uint64_t*)myalloc(((lenA+lenB+w-1+w)/w)*sizeof(uint64_t));
     size = (lenA + lenB)/4 + 1;
     Q = (queue2*)myalloc(size*sizeof(queue2));
     ptr = ptrA = ptrB = 0;
     elems = 0;
     head = 0; tail = 1; prev = 0;
     Q[0] = hasdist(1|2,1,level-1);
     if (levels) for (v=0;v<level;v++) levels[v] = 0;
     while (head != tail)
        { uint has = gethas(Q[head]);
	  uint dist = getdist(Q[head]);
   	  uint num,ones;
	  if (has == 1) // copy node of A
	     { num = 4 * getnumc(Q[head]);
	       if (levels) levels[dist] += num/4;
	       ones = 0;
	       while (num)
		 { n = mmin(num,w);
		   num -= n;
	           sigs = bitsReadA(treeA,ptrA,n);
	           ptrA += n;
	           ones += popcount(sigs);
	           addWord(treeM,ptr,sigs);
	           ptr += n;
		 }
	       if (dist) // don't add tasks for leaves
	          { while (ones)
		      { o = mmin(ones,256); // it is stored in a byte
		        ones -= o;
			if ((prev != head) && (gethas(Q[prev]) == has) &&
			    (getdist(Q[prev]) == dist-1) &&
			    (getnumc(Q[prev])+o <= 256))
			       addnumc(Q[prev],o);
		        else { Q[tail] = hasdist(has,o,dist-1);
		               prev = tail; tail = (tail+1)%size;
			     }
		      }
		  }
	       else elems += ones;
	     }
	  else if (has == 2) // copy node of B
	     { num = 4 * getnumc(Q[head]);
	       if (levels) levels[dist] += num/4;
	       ones = 0;
	       while (num)
		 { n = mmin(num,w);
		   num -= n;
	           sigs = bitsReadA(treeB,ptrB,n);
	           ptrB += n;
	           ones += popcount(sigs);
	           addWord(treeM,ptr,sigs);
	           ptr += n;
		 }
	       if (dist) // don't add tasks for leaves
	          { while (ones)
		      { o = mmin(ones,w/4);
		        ones -= o;
			if ((prev != head) && (gethas(Q[prev]) == has) &&
			    (getdist(Q[prev]) == dist-1) &&
			    (getnumc(Q[prev])+o <= w/4))
			       addnumc(Q[prev],o);
		        else { Q[tail] = hasdist(has,o,dist-1);
		               prev = tail; tail = (tail+1)%size;
			     }
		      }
		  }
	       else elems += ones;
	     }
	  else // both A and B nodes exist
	     { sigA = (treeA[ptrA/w] >> (ptrA%w)) & 0xF;
	       sigB = (treeB[ptrB/w] >> (ptrB%w)) & 0xF;
	       ptrA += 4; ptrB += 4;
	       sig = sigA | sigB;
	       if (dist)
	          for (v=0;v<4;v++)
	              { uint po;
			if (sig & (1<<v))
		           { has = ((sigA >> v) & 1)*1 + ((sigB >> v) & 1)*2;
			     if ((has != (1|2)) && (prev != head) &&
				 (gethas(Q[prev]) == has) &&
			         (getdist(Q[prev]) == dist-1) &&
			         (getnumc(Q[prev])+1 <= w/4))
			          addnumc(Q[prev],1);
			     else
		                { Q[tail] = hasdist(has,1,dist-1);
		                  prev = tail; tail = (tail+1)%size;
		                }
			   }
		       }
	       else elems += pop4(sig);
	       treeM[ptr/w] &= ~(((uint64_t)0xF) << (ptr%w));
	       treeM[ptr/w] |= ((uint64_t)sig) << (ptr%w);
	       ptr += 4;
	       if (levels) levels[dist]++;
	     }
	  head = (head+1)%size;
        }
     myfree(Q);
     treeM = (uint64_t*)myrealloc(treeM,((ptr+w-1)/w)*sizeof(uint64_t));
     *telems = elems;
     *len = ptr;
     return treeM;
   }

