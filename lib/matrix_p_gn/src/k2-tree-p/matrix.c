
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "k2-tree-p/utilstime.h"
#include "k2-tree-p/matrix.h"

static uint *mapA,*mapB;

typedef struct {
   uint64_t *tree,*levels;
   uint64_t elems,len;
   } partition;

static partition empty = { NULL, NULL, 0, 0 };

uint DLEV; // only first DLEV levels of matrices are parallel
uint OLEV; // last OLEV levels of matrices are sequential

#define PAR(level,dlevel) (((dlevel) < DLEV) && ((level) >= OLEV))

typedef struct
   { k2tree treeA,treeB;
     k2node nodeA,nodeB;
     k2node nodeA1,nodeB1,nodeA2,nodeB2;
     uint level, dlevel, plevel;
     partition *res;
     k2tree *tres;
     uint64_t row,col;
   } Tpar;

	// creates matrix of width x height with n cells (2n ints row,col) 
	// reorders cells array

matrix matCreate (uint64_t height, uint64_t width, uint64_t n, uint64_t *cells)

   { matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(mmax(width,height)-1);
     mat->elems = n;
     mat->transposed = 0;
     if (n == 0) mat->tree = NULL;
     else mat->tree = k2create(n,mat->logside,cells);
     return mat;
   }
     
matrix matCreate32 (uint64_t height, uint64_t width, uint64_t n, uint32_t *cells)

   { matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(mmax(width,height)-1);
     mat->elems = n;
     mat->transposed = 0;
     if (n == 0) mat->tree = NULL;
     else mat->tree = k2create32(n,mat->logside,cells);
     return mat;
   }
     
        // destroys matrix

void matDestroy (matrix M)

   { if (M == NULL) return; // safety
     k2destroy (M->tree);
     myfree(M);
   }

        // creates an empty matrix

matrix matEmpty (uint64_t height, uint64_t width)

   { matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(mmax(width,height)-1);
     mat->elems = 0;
     mat->transposed = 0;
     mat->tree = NULL;
     return mat;
   }

        // creates a matrix with cell row,col

matrix matOne (uint64_t height, uint64_t width, uint64_t row, uint64_t col)

   { matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
     uint64_t *data;
     uint64_t sig,ptr;
     int bit,len,words;
     mat->width = width;
     mat->height = height;
     len = mat->logside = numbits(mmax(width,height)-1);
     mat->elems = 1;
     mat->transposed = 0;
     words = (4*len+w-1)/w;
     data = (uint64_t*)myalloc(words*sizeof(uint64_t));
     for (ptr=0;ptr<words;ptr++) data[ptr] = 0;
     ptr = 0;
     for (bit=len-1;bit>=0;bit--)
	 { sig = 1 << (2*((row >> bit) & 1) + ((col >> bit) & 1));
	   data[ptr/w] |= sig << (ptr%w);
	   ptr += 4;
	 }
     mat->tree = k2createFrom(len,ptr,data,1);
     return mat;
   }

        // creates an identity matrix

matrix matId (uint64_t side)

   { matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
     uint64_t *data;
     uint64_t p,k;
     uint *levels,*odd;
     uint level,l;

     side = side ? side : 1;
     mat->width = side;
     mat->height = side;
     mat->logside = level = numbits(side-1);
     mat->elems = side;
     mat->transposed = 0;
     p = (((uint64_t)1) << level)-1; // upper bound to allocate
     data = (uint64_t*)myalloc(((4*p+w-1)/w)*sizeof(uint64_t));
     memset(data,0x99,((4*p+w-1)/w)*sizeof(uint64_t));
     levels = (uint*)myalloc((1+level)*sizeof(uint));
     odd = (uint*)myalloc((1+level)*sizeof(uint));
     k = side;
     for (l=1;l<=level;l++)
         { levels[l] = (k+1)/2; odd[l] = k % 2;
           k = (k+1)/2;
         }
     p = 0;
     for (l=level;l>=1;l--)
         { p += levels[l];
           if (odd[l]) bitsWriteA(data,4*(p-1)+3,0);
         }
     myfree(levels); myfree(odd);
     data = (uint64_t*)myrealloc(data,((4*p+w-1)/w)*sizeof(uint64_t));
     mat->tree = k2createFrom(level,4*p,data,1);
     return mat;
   }

	// creates a new copy of A, with its own data

matrix matCopy (matrix A)

   { matrix M = (matrix)myalloc(sizeof(struct s_matrix));
     *M = *A;
     if (A->tree != NULL) M->tree = k2copy(A->tree);
     return M;
   }

        // transpose a matrix, creating a non-allocated copy that shares the 
        // data. You need not (and should not) matDestroy this copy
	// use *M = matTranspose(M) to actually transpose M
struct s_matrix matTranspose (matrix M)

   { struct s_matrix T = *M;
     T.transposed = 1-T.transposed;
     return T;
   }

        // writes M to file, which must be opened for writing

void matSave (matrix M, FILE *file)

   { uint aux = M->logside + (M->transposed << 16);;
     fwrite (&M->elems,sizeof(uint64_t),1,file);
     fwrite (&aux,sizeof(uint),1,file);
     fwrite (&M->width,sizeof(uint64_t),1,file);
     fwrite (&M->height,sizeof(uint64_t),1,file);
     if (M->elems != 0) k2save (M->tree,file);
   }

        // loads matrix from file, which must be opened for reading

matrix matLoad (FILE *file)

    { matrix M;
      uint aux;
      M = (matrix)myalloc(sizeof(struct s_matrix));
      fread (&M->elems,sizeof(uint64_t),1,file);
      fread (&aux,sizeof(uint),1,file);
      M->logside = aux & ((1<<16)-1);
      M->transposed = aux >> 16;
      fread (&M->width,sizeof(uint64_t),1,file);
      fread (&M->height,sizeof(uint64_t),1,file);
      if (M->elems == 0) M->tree = NULL;
      else M->tree = k2load(file);
      return M;
    }

	// space of the matrix, in w-bit words

uint64_t matSpace (matrix M)

   { uint64_t space = sizeof(struct s_matrix)*8/w;
     if (M == NULL) return 0; // safety
     if (M->tree != NULL) space += k2space(M->tree);
     return space;
   }

        // dimensions of M, returns elems

uint64_t matDims (matrix M, uint *logside, uint64_t *width, uint64_t *height)

   { if (logside != NULL) *logside = M->logside;
     if (width != NULL) *width = M->transposed ? M->height: M->width;
     if (height != NULL) *height = M->transposed ? M->width : M->height;
     return M->elems;
   }

        // accesses a cell

uint matAccess (matrix M, uint64_t row, uint64_t col)

   { k2tree T = M->tree;
     k2node node;
     uint level;
     uint64_t aux;
     if (M->elems == 0) return 0;
     node = k2root(T);
     level = k2levels(T)-1;
     if (M->transposed) { aux = row; row = col; col = aux; }
     while (1)
       { uint v = 2 * ((row >> level) & 1) + ((col >> level) & 1);
	 if (!k2hasChild(T,node,v)) return 0;
	 if (level == 0) return 1;
	 node = k2child(T,node,v);
	 level--;
       }
   }

        // recovers all the cells in [r1..r2] x [c1..c2]
        // writes 2n integers in buffer, which must be preallocated
        // to 2*elems of uint64_t (worst case)
	// returns number of cells. just counts if buffer is NULL

uint64_t matCollect (matrix M, uint64_t r1, uint64_t r2, 
		     uint64_t c1, uint64_t c2, uint64_t *buffer)

   { uint64_t aux;
     if (M->elems == 0) return 0;
     if (M->transposed)
	{ aux = c1; c1 = r1; r1 = aux;
	  aux = c2; c2 = r2; r2 = aux;
	}
     if (r2 == fullSide) r2 = M->height-1; // already in concrete repres
     if (c2 == fullSide) c2 = M->width-1;
     return k2collect (M->tree,r1,r2,c1,c2,buffer,M->transposed);
   }

uint64_t matCollect32 (matrix M, uint32_t r1, uint32_t r2, 
		     uint32_t c1, uint32_t c2, uint32_t *buffer)

   { uint64_t aux;
     if (M->elems == 0) return 0;
     if (M->transposed)
	{ aux = c1; c1 = r1; r1 = aux;
	  aux = c2; c2 = r2; r2 = aux;
	}
     if (r2 == fullSide32) r2 = M->height-1; // already in concrete repres
     if (c2 == fullSide32) c2 = M->width-1;
     return k2collect32 (M->tree,r1,r2,c1,c2,buffer,M->transposed);
   }

	// (boolean) sum of two matrices, assumed to be of the same logside

matrix matSum (matrix A, matrix B)

   { uint64_t *sum;
     uint64_t len,elems;
     uint64_t *treeA,*treeB;
     uint64_t width,height;
     uint64_t lenA,lenB;
     matrix M;

     if (A->logside != B->logside)
	{ fprintf(stderr,"Error: sum of matrices of different side\n");
	  exit(1);
	}
     if (A->elems == 0) M = matCopy(B);
     else if (B->elems == 0) M = matCopy(A);
     else { M = (matrix)myalloc(sizeof(struct s_matrix));
            M->logside = A->logside;
            treeA = bitsData(k2bits(A->tree)); 
	    treeB = bitsData(k2bits(B->tree));
            lenA = bitsLength(k2bits(A->tree)); 
	    lenB = bitsLength(k2bits(B->tree));
	    if (A->transposed == B->transposed) // transposition can be ignored
               { sum = k2merge (treeA,lenA,treeB,lenB,k2levels(A->tree),
		                &len,&elems,NULL);
		 M->transposed = A->transposed;
	       }
	    else if (A->elems > B->elems) // better that B is transposed
	       { sum = k2mergeT(treeA,lenA,B->tree,&len,&elems,NULL);
		 M->transposed = A->transposed;
	       }
	    else  // better that A is transposed
	       { sum = k2mergeT(treeB,lenB,A->tree,&len,&elems,NULL);
		 M->transposed = B->transposed;
	       }
            M->elems = elems;
            M->tree = k2createFrom (k2levels(A->tree),len,sum,1);
	  }
     if (A->transposed == B->transposed)
	{ M->width = mmax(A->width,B->width);
     	  M->height = mmax(A->height,B->height);
	}
     else if (A->elems > B->elems) // the transposition of A dominates
	{ M->width = mmax(A->width,B->height);
     	  M->height = mmax(A->height,B->width);
	}
     else // A->elems <= B->elems, the transposition of B dominates
	{ M->width = mmax(B->width,A->height);
     	  M->height = mmax(B->height,A->width);
	}
     return M;
   }

        // version with one row or one column, or both

	// copies len bits starting at *src + psrc 
	// to tgt from bit position ptgt 
	// WARNING: leave at least one extra word to spare in tgt

static void copyBits (uint64_t *tgt, uint64_t ptgt, 
		      uint64_t *src, uint64_t psrc, uint64_t len)

   { uint64_t old,mask,shift;

	// easy cases if they are similarly misaligned
	// I didn't align to 8 due to possible endianness issues
     tgt += ptgt/w; ptgt %= w;
     src += psrc/w; psrc %= w;
     if (ptgt == psrc) 
	{ if (ptgt != 0)
	     { shift = w-ptgt;
	       mask = (~((uint64_t)0)) >> shift;
               *tgt = (*tgt & mask) + (*src & ~mask);
	       *tgt++; *src++; len -= shift;
	     }
	  memcpy (tgt,src,((len+w-1)/w)*sizeof(uint64_t)); 
	  return;
	}
	// general case, we first align the source
     mask = (((uint64_t)1) << ptgt)-1;
     if (ptgt < psrc) // first word from src fits in ptgt
	{ *tgt = (*tgt & mask) | ((*src++ >> (psrc-ptgt)) & ~mask);
	  ptgt += w-psrc;
	}
     else
	{ *tgt = (*tgt & mask) | ((*src << (ptgt-psrc)) & ~mask);
	  ptgt -= psrc;
	  *++tgt = *src++ >> (w-ptgt);
	}
     if (len <= w-psrc) return;
     len -= w-psrc;
     	// now src is aligned, copy all the rest
     mask = (((uint64_t)1) << ptgt)-1;
     old = *tgt & mask;
     len += w;
     while (len >= w)
	{ *tgt++ = old | (*src << ptgt);
	  old = *src++ >> (w-ptgt);
	  len -= w;
	}
     *tgt = old;
   }

static partition compose (partition *part, uint level, uint freepart)

   { uint64_t ptr,old;
     uint l,v,sig;
     partition answer;
     uint64_t ptrs[4];

     answer.elems = part[0].elems+part[1].elems+part[2].elems+part[3].elems;
     if (answer.elems == 0) return empty;
     sig =  (part[0].elems != 0) +
           ((part[1].elems != 0) << 1) +
           ((part[2].elems != 0) << 2) +
           ((part[3].elems != 0) << 3);
     if (pop4(sig) == 1) // trivial concatenation
        { v = (((sig >> 1) & 1) * 1) +
              (((sig >> 2) & 1) * 2) +
              (((sig >> 3) & 1) * 3);
          answer.len = 4 + part[v].len;
          answer.levels = (uint64_t*)myrealloc(part[v].levels,
                                             level*sizeof(uint64_t));
          answer.levels[level-1] = 1;
          answer.tree = (uint64_t*)myalloc((1+(answer.len+w-1)/w)
                                           *sizeof(uint64_t)); //+1 for copyBits
          answer.tree[0] = sig;
          copyBits(answer.tree,4,part[v].tree,0,part[v].len);
          if (freepart) myfree(part[v].tree);
	  part[v].levels = NULL;
          return answer;
        }
        // nontrivial concatenation
     answer.len = 4 + part[0].len+part[1].len+part[2].len+part[3].len;
     answer.levels = (uint64_t*)myalloc(level*sizeof(uint64_t));
     answer.levels[level-1] = 1;
     answer.tree = (uint64_t*)myalloc((1+(answer.len+w-1)/w)
                                         *sizeof(uint64_t)); //+1 for copyBits
     answer.tree[0] = sig;
     for (v=0;v<4;v++) ptrs[v] = 0;
     ptr = old = 4;
     for (l=level-1;l>0;)
         { uint64_t segm;
           l--;
           for (v=0;v<4;v++)
               { if (part[v].tree != NULL)
                    { segm = 4*part[v].levels[l];
                      copyBits(answer.tree,ptr,part[v].tree,ptrs[v],segm);
                      ptr += segm; ptrs[v] += segm;
                    }
               }
           answer.levels[l] = (ptr - old)/4; old = ptr;
         }
     if (freepart)
        for (v=0;v<4;v++)
            { myfree(part[v].tree);
              myfree(part[v].levels);
            }
     return answer;
   }

#define valA(i) ((sigA >> i) & 1)
#define valB(i) ((sigB >> i) & 1)
#define val(i) ((sig >> i) & 1)

static partition single (uint sig)

   { partition answer;
     if (sig == 0) return empty;
     answer.elems = pop4(sig);
     answer.len = 4;
     answer.tree = (uint64_t*)myalloc(((4+w-1)/w)*sizeof(uint64_t));
     answer.tree[0] = sig;
     answer.levels = (uint64_t*)myalloc(1*sizeof(uint64_t));
     answer.levels[0] = 1;
     return answer;
   }

static inline void fixSig (uint64_t row, uint64_t col, uint64_t lim,
                      uint *sigA, uint *sigB)

   { if (row != fullSide)
        { if (row >= lim) { *sigA &= ~0x3; *sigB &= ~0x3; }
          else { *sigA &= ~0xC; *sigB &= ~0xC; }
        }
     if (col != fullSide)
        { if (col >= lim) { *sigA &= ~0x5; *sigB &= ~0x5; }
          else { *sigA &= ~0xA; *sigB &= ~0xA; }
        }
   }

static partition k2sumRC (k2tree treeA,k2node nodeA,k2tree treeB,k2node nodeB,
			  uint level, uint64_t row, uint64_t col, uint has)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     uint64_t lim;

     if (has & 1) sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     else sigA = 0;
     if (has & 2) sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);
     else sigB = 0;

     lim = ((uint64_t)1)<<(level-1);
     fixSig(row,col,lim,&sigA,&sigB);

        // base case, level = 1
     if (level == 1) return single (sigA | sigB);

     for (v=0;v<4;v++) part[v] = empty;
     has = valA(0) + 2*valB(0);
     if (has) part[0] = k2sumRC (treeA,childA[0],treeB,childB[0],
                                 level-1,row,col,has);
     has = valA(1) + 2*valB(1);
     if (has) part[1] = k2sumRC (treeA,childA[1],treeB,childB[1],
                              level-1,row,col == fullSide ? col : col-lim, has);
     has = valA(2) + 2*valB(2);
     if (has) part[2] = k2sumRC (treeA,childA[2],treeB,childB[2],
                              level-1,row == fullSide ? row : row-lim,col, has);
     has = valA(3) + 2*valB(3);
     if (has) part[3] = k2sumRC (treeA,childA[3],treeB,childB[3],
                                 level-1,row == fullSide ? row : row-lim,
                                         col == fullSide ? col : col-lim, has);

        // combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static uint64_t *k2sum1 (k2tree A, k2tree B, uint64_t row, uint64_t col,
			 uint64_t *len, uint64_t *telems)

   { partition sum;
     uint level,has;

     if ((A == NULL) && (B == NULL))
	{ *len = *telems = 0; return NULL; }
     has = 0;
     if (A != NULL) { level = k2levels(A); has += 1; }
     if (B != NULL) { level = k2levels(B); has += 2; }
     
     sum = k2sumRC(A,k2root(A),B,k2root(B),level,row,col,has);

     *telems = sum.elems;
     *len = sum.len;
     myfree (sum.levels);
     return sum.tree;
   }

matrix matSum1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { uint64_t *sum;
     matrix M;
     uint64_t wA,hA,wB,hB;
     uint64_t len;

     if ((row == fullSide) && (col == fullSide)) return matSum(A,B);

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: sum of matrices of different side\n");
          exit(1);
        }
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);

     if ((row != fullSide) && (col != fullSide)) // just one cell
        { if (matAccess(A,row,col) || matAccess(B,row,col))
	     return matOne(mmax(hA,hB),mmax(wA,wB),row,col);
	  else return matEmpty(mmax(hA,hB),mmax(wA,wB));
        }

     M = (matrix)myalloc(sizeof(struct s_matrix));
     M->logside = A->logside;
     M->width = mmax(wA,wB); M->height = mmax(hA,hB);
     M->transposed = 0;
     mapA = A->transposed ? mapTr : mapId;
     mapB = B->transposed ? mapTr : mapId;
     sum = k2sum1 (A->tree,B->tree,row,col,&len,&M->elems);
     if (M->elems == 0) M->tree = NULL;
     else if (A->tree) M->tree = k2createFrom (k2levels(A->tree),len,sum,1);
     else M->tree = k2createFrom (k2levels(B->tree),len,sum,1);
     return M;
   }

	// (boolean) difference A - B, assumed to be of the same logside

	// extracts the subtree as a partition, allows transpositions

static partition k2extract (k2tree treeA, k2node nodeA, uint level)

   { partition part[4];
     k2node childA[4];
     uint sigA;
     uint v;

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);

        // base case, level = 1
     if (level == 1) return single (sigA);

     for (v=0;v<4;v++)
         { if (valA(v)) part[v] = k2extract(treeA,childA[v],level-1);
           else part[v] = empty;
         }
        // combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parDif (Tpar *m);

static partition k2dif (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		        uint level, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

        // base case, level = 1
     if (level == 1) return single (sigA & ~sigB);

     	// compute the 4 quadrants
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
	     { if (!valB(v)) part[v] = k2extract(treeA,childA[v],level-1);
               else part[v] = k2dif(treeA,childA[v],treeB,childB[v],
			            level-1,dlevel+1);
	     }
	  else
	     { m[v].treeA = treeA; m[v].nodeA = childA[v];
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parDif,&m[v]))
		  { parDif(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parDif (Tpar *m)
   { if (m->treeB == NULL)
	  *(m->res) = k2extract(m->treeA,m->nodeA,m->level);
     else *(m->res) = k2dif(m->treeA,m->nodeA,m->treeB,m->nodeB,
		     	    m->level,m->dlevel);
   }

static void parDifRC (Tpar *m);

static partition k2difRC (k2tree treeA, k2node nodeA, k2tree treeB,k2node nodeB,
		          uint level, uint64_t row, uint64_t col, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     uint64_t lim;
     uint64_t rows[4],cols[4];
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

     lim = ((uint64_t)1)<<(level-1);
     fixSig(row,col,lim,&sigA,&sigB);

        // base case, level = 1
     if (level == 1) return single (sigA & ~sigB);

     rows[0] = rows[1] = row; cols[0] = cols[2] = col;
     rows[2] = rows[3] = row == fullSide ? row : row-lim;
     cols[1] = cols[3] = col == fullSide ? col : col-lim;
	// compute the 4 quadrants
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
	     { if (!valB(v)) part[v] = k2sumRC(treeA,childA[v],NULL,0,
			  		       level-1,rows[v],cols[v],1);
               else part[v] = k2difRC(treeA,childA[v],treeB,childB[v],
			  	 level-1,rows[v],cols[v],dlevel+1);
	     }
	  else
	     { m[v].treeA = treeA; m[v].nodeA = childA[v];
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].row = rows[v]; m[v].col = cols[v];
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parDifRC,&m[v]))
		  { parDifRC(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parDifRC (Tpar *m)
   { if (m->treeB == NULL)
	  *(m->res) = k2sumRC(m->treeA,m->nodeA,NULL,0,
		              m->level,m->row,m->col,1);
     else *(m->res) = k2difRC(m->treeA,m->nodeA,m->treeB,m->nodeB,
		     	      m->level,m->row,m->col,m->dlevel);
   }

matrix matDif (matrix A, matrix B)

   { return matDif1 (fullSide,A,B,fullSide);
   }

matrix matDif1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { partition dif;
     matrix M;
     uint64_t wA,hA,wB,hB,wM,hM;

     DLEV = 7; OLEV = 5;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: diff of matrices of different logside\n");
          exit(1);
        }
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
     wM = mmax(wA,wB); hM = mmax(hA,hB);
     if (A->elems == 0) return matEmpty(hM,wM);
     if (B->elems == 0)
	{ M = matSum1(row,A,B,col); // has provisions for this case
	}
     else
        { mapA = A->transposed ? mapTr : mapId;
          mapB = B->transposed ? mapTr : mapId;
          M = (matrix)myalloc(sizeof(struct s_matrix));
          M->logside = A->logside;
	  if ((row == fullSide) && (col == fullSide))
               dif = k2dif (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		   	    k2levels(A->tree),0);
	  else dif = k2difRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		              k2levels(A->tree),row,col,0);
          myfree (dif.levels);
          M->elems = dif.elems;
          if (M->elems == 0) M->tree = NULL;
          else M->tree = k2createFrom (k2levels(A->tree),dif.len,dif.tree,1);
          M->transposed = 0;
	}
     if (M->transposed) { M->height = wM; M->width = hM; }
     else { M->height = hM; M->width = wM; }
     return M;
   }

	// (boolean) symmetric difference, assumed to be of the same logside

static void parXor (Tpar *m);

static partition k2xor (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		        uint level, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

        // base case, level = 1
     if (level == 1) return single (sigA ^ sigB);

     	// compute the 4 quadrants
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) && !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
	     { if (!valB(v)) part[v] = k2extract(treeA,childA[v],level-1);
	       else if (!valA(v)) part[v] = k2extract(treeB,childB[v],level-1);
               else part[v] = k2xor(treeA,childA[v],treeB,childB[v],
			            level-1,dlevel+1);
	     }
	  else
	     { if (valA(v)) { m[v].treeA = treeA; m[v].nodeA = childA[v]; }
	       else m[v].treeA = NULL;
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parXor,&m[v]))
		  { parXor(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parXor (Tpar *m)
   { if (m->treeB == NULL)
	  *(m->res) = k2extract(m->treeA,m->nodeA,m->level);
     else if (m->treeA == NULL)
	  *(m->res) = k2extract(m->treeB,m->nodeB,m->level);
     else *(m->res) = k2xor(m->treeA,m->nodeA,m->treeB,m->nodeB,
		     	    m->level,m->dlevel);
   }

static void parXorRC (Tpar *m);

static partition k2xorRC(k2tree treeA, k2node nodeA, k2tree treeB,k2node nodeB,
		         uint level, uint64_t row, uint64_t col, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     uint64_t lim;
     uint64_t rows[4],cols[4];
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

     lim = ((uint64_t)1)<<(level-1);
     fixSig(row,col,lim,&sigA,&sigB);

        // base case, level = 1
     if (level == 1) return single (sigA ^ sigB);

     rows[0] = rows[1] = row; cols[0] = cols[2] = col;
     rows[2] = rows[3] = row == fullSide ? row : row-lim;
     cols[1] = cols[3] = col == fullSide ? col : col-lim;
	// compute the 4 quadrants
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) && !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
	     { if (!valB(v)) part[v] = k2sumRC(treeA,childA[v],NULL,0,
			  		       level-1,rows[v],cols[v],1);
	       else if (!valA(v)) part[v] = k2sumRC(treeB,childB[v],NULL,0,
			  		       level-1,rows[v],cols[v],1);
               else part[v] = k2xorRC(treeA,childA[v],treeB,childB[v],
			  	      level-1,rows[v],cols[v],dlevel+1);
	     }
	  else
	     { if (valA(v)) { m[v].treeA = treeA; m[v].nodeA = childA[v]; }
	       else m[v].treeA = NULL;
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].row = rows[v]; m[v].col = cols[v];
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parXorRC,&m[v]))
		  { parXorRC(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parXorRC (Tpar *m)
   { if (m->treeB == NULL)
	  *(m->res) = k2sumRC(m->treeA,m->nodeA,NULL,0,
		              m->level,m->row,m->col,1);
     else if (m->treeA == NULL)
	  *(m->res) = k2sumRC(m->treeB,m->nodeB,NULL,0,
		              m->level,m->row,m->col,1);
     else *(m->res) = k2xorRC(m->treeA,m->nodeA,m->treeB,m->nodeB,
		     	      m->level,m->row,m->col,m->dlevel);
   }

matrix matXor (matrix A, matrix B)

   { return matXor1 (fullSide,A,B,fullSide);
   }

matrix matXor1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { partition xor;
     matrix M;
     uint64_t wA,hA,wB,hB,wM,hM;

     DLEV = 7; OLEV = 5;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: sym diff of matrices of different logside\n");
          exit(1);
        }
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
     wM = mmax(wA,wB); hM = mmax(hA,hB);
     if ((A->elems == 0) && (B->elems == 0)) return matEmpty(hM,wM);
     if ((A->elems == 0) || (B->elems == 0))
        { M = matSum1(row,A,B,col); // has provisions for this case
        }
     else
        { mapA = A->transposed ? mapTr : mapId;
          mapB = B->transposed ? mapTr : mapId;
          M = (matrix)myalloc(sizeof(struct s_matrix));
          M->logside = A->logside;
	  if ((row == fullSide) && (col == fullSide))
               xor = k2xor (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		   	    k2levels(A->tree),0);
	  else xor = k2xorRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		              k2levels(A->tree),row,col,0);
          myfree (xor.levels);
          M->elems = xor.elems;
          if (M->elems == 0) M->tree = NULL;
          else M->tree = k2createFrom (k2levels(A->tree),xor.len,xor.tree,1);
          M->transposed = 0;
	}
     if (M->transposed) { M->height = wM; M->width = hM; }
     else { M->height = hM; M->width = wM; }
     return M;
   }

	// (boolean) intersection of A and B, of the same size

static void parAnd (Tpar *m);

static partition k2and (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		        uint level, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

        // base case, level = 1
     if (level == 1) return single (sigA & sigB);

     	// compute the 4 quadrants, could be parallelized
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) || !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
             { part[v]= k2and(treeA,childA[v],treeB,childB[v],level-1,dlevel+1);
	     }
	  else
	     { m[v].treeA = treeA; m[v].nodeA = childA[v];
	       m[v].treeB = treeB; m[v].nodeB = childB[v];
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parAnd,&m[v]))
		  { parAnd(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parAnd (Tpar *m)
   { *(m->res) = k2and(m->treeA,m->nodeA,m->treeB,m->nodeB,m->level,m->dlevel);
   }

static void parAndRC (Tpar *m);

static partition k2andRC (k2tree treeA, k2node nodeA, k2tree treeB,k2node nodeB,
		          uint level, uint64_t row, uint64_t col, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     uint64_t lim;
     uint64_t rows[4],cols[4];
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

     lim = ((uint64_t)1)<<(level-1);
     fixSig(row,col,lim,&sigA,&sigB);

        // base case, level = 1
     if (level == 1) return single (sigA & sigB);

     rows[0] = rows[1] = row; cols[0] = cols[2] = col;
     rows[2] = rows[3] = row == fullSide ? row : row-lim;
     cols[1] = cols[3] = col == fullSide ? col : col-lim;
	// compute the 4 quadrants (could be parallelized)
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) || !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
             { part[v] = k2andRC(treeA,childA[v],treeB,childB[v],
			  	 level-1,rows[v],cols[v],dlevel+1);
	     }
	  else
	     { m[v].treeA = treeA; m[v].nodeA = childA[v];
	       m[v].treeB = treeB; m[v].nodeB = childB[v];
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].row = rows[v]; m[v].col = cols[v];
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parAndRC,&m[v]))
		  { parAndRC(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parAndRC (Tpar *m)
   { *(m->res) = k2andRC (m->treeA,m->nodeA,m->treeB,m->nodeB,
		          m->level,m->row,m->col,m->dlevel);
   }

matrix matAnd (matrix A, matrix B)

   { return matAnd1 (fullSide,A,B,fullSide);
   }

matrix matAnd1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { partition and;
     matrix M;
     uint64_t wA,hA,wB,hB,wM,hM;

     //DLEV = 7; OLEV = 5;
     DLEV = 2; OLEV = 1;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: and of matrices of different logside\n");
          exit(1);
        }
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
     wM = mmax(wA,wB); hM = mmax(hA,hB);
     if ((A->elems == 0) || (B->elems == 0)) return matEmpty(hM,wM);
     M = (matrix)myalloc(sizeof(struct s_matrix));
     M->logside = A->logside;
     mapA = A->transposed ? mapTr : mapId;
     mapB = B->transposed ? mapTr : mapId;
     if ((row == fullSide) && (col == fullSide))
          and = k2and (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		     k2levels(A->tree),0);
     else and = k2andRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		         k2levels(A->tree),row,col,0);
     myfree (and.levels);
     M->elems = and.elems;
     if (M->elems == 0) M->tree = NULL;
     else M->tree = k2createFrom (k2levels(A->tree),and.len,and.tree,1);
     M->transposed = 0;
     M->height = hM; M->width = wM;
     return M;
   }

	// (boolean) union of A and B, of the same size
	// same result as matSum but this relies on parallelism to be faster

static void parOr (Tpar *m);

static partition k2or (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		       uint level, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

        // base case, level = 1
     if (level == 1) return single (sigA | sigB);

     	// compute the 4 quadrants, could be parallelized
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) && !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
	     { if (!valA(v)) part[v] = k2extract(treeB,childB[v],level-1);
	       else if (!valB(v)) part[v] = k2extract(treeA,childA[v],level-1);
	       else part[v] = k2or(treeA,childA[v],treeB,childB[v],
			           level-1,dlevel+1);
	     }
	  else
	     { if (valA(v)) { m[v].treeA = treeA; m[v].nodeA = childA[v]; }
	       else m[v].treeA = NULL;
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parOr,&m[v]))
		  { parOr(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parOr (Tpar *m)
   { if (m->treeA == NULL)
	  *(m->res) = k2extract(m->treeB,m->nodeB,m->level);
     else if (m->treeB == NULL)
	  *(m->res) = k2extract(m->treeA,m->nodeA,m->level);
     else *(m->res) = k2or(m->treeA,m->nodeA,m->treeB,m->nodeB,
		           m->level,m->dlevel);
   }

static void parOrRC (Tpar *m);

static partition k2orRC (k2tree treeA, k2node nodeA, k2tree treeB,k2node nodeB,
		         uint level, uint64_t row, uint64_t col, uint dlevel)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;
     uint64_t lim;
     uint64_t rows[4],cols[4];
     pthread_t thr[4];
     Tpar m[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

     lim = ((uint64_t)1)<<(level-1);
     fixSig(row,col,lim,&sigA,&sigB);

        // base case, level = 1
     if (level == 1) return single (sigA | sigB);

     rows[0] = rows[1] = row; cols[0] = cols[2] = col;
     rows[2] = rows[3] = row == fullSide ? row : row-lim;
     cols[1] = cols[3] = col == fullSide ? col : col-lim;
	// compute the 4 quadrants (could be parallelized)
     for (v=0;v<4;v++)
        { thr[v] = 0;
	  if (!valA(v) && !valB(v)) part[v] = empty;
          else if (!PAR(level,dlevel))
             { if (!valA(v)) part[v] = k2sumRC(treeB,childB[v],NULL,0,
			     		       level-1,rows[v],cols[v],1);
	     else if (!valB(v)) part[v] = k2sumRC(treeA,childA[v],NULL,0,
			     		          level-1,rows[v],cols[v],1);
	     else part[v] = k2orRC(treeA,childA[v],treeB,childB[v],
			  	   level-1,rows[v],cols[v],dlevel+1);
	     }
	  else
	     { if (valA(v)) { m[v].treeA = treeA; m[v].nodeA = childA[v]; }
	       else m[v].treeA = NULL;
	       if (valB(v)) { m[v].treeB = treeB; m[v].nodeB = childB[v]; }
	       else m[v].treeB = NULL;
	       m[v].level = level-1; m[v].dlevel = dlevel+1;
	       m[v].row = rows[v]; m[v].col = cols[v];
	       m[v].res = &part[v];
	       if (pthread_create(&thr[v],NULL,(void*(*)(void*))parOrRC,&m[v]))
		  { parOrRC(&m[v]); thr[v] = 0; }
	     }
	}
     for (v=0;v<4;v++)
	 if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parOrRC (Tpar *m)
   { if (m->treeA == NULL)
	  *(m->res) = k2sumRC(m->treeB,m->nodeB,NULL,0,
			      m->level,m->row,m->col,1);
     else if (m->treeB == NULL)
	  *(m->res) = k2sumRC(m->treeA,m->nodeA,NULL,0,
			      m->level,m->row,m->col,1);
     else *(m->res) = k2orRC(m->treeA,m->nodeA,m->treeB,m->nodeB,
		             m->level,m->row,m->col,m->dlevel);
   }


matrix matOr (matrix A, matrix B)

   { return matOr1 (fullSide,A,B,fullSide);
   }

matrix matOr1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { partition or;
     matrix M;
     uint64_t wA,hA,wB,hB,wM,hM;

     //DLEV = 6; OLEV = 5;
     DLEV = 2; OLEV = 1;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: or of matrices of different logside\n");
          exit(1);
        }
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
     wM = mmax(wA,wB); hM = mmax(hA,hB);
     if ((A->elems == 0) && (B->elems == 0)) return matEmpty(hM,wM);
     if ((A->elems == 0) || (B->elems == 0))
        { M = matSum1(row,A,B,col); // has provisions for this case
        }
     else
	{ M = (matrix)myalloc(sizeof(struct s_matrix));
          M->logside = A->logside;
          mapA = A->transposed ? mapTr : mapId;
          mapB = B->transposed ? mapTr : mapId;
          if ((row == fullSide) && (col == fullSide))
               or = k2or (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		          k2levels(A->tree),0);
          else or = k2orRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		            k2levels(A->tree),row,col,0);
          myfree (or.levels);
          M->elems = or.elems;
          M->tree = k2createFrom (k2levels(A->tree),or.len,or.tree,1);
          M->transposed = 0;
	}
     if (M->transposed) { M->height = wM; M->width = hM; }
     else { M->height = hM; M->width = wM; }
     return M;
   }


	// (boolean) product of two matrices, assumed to be of the same side

static void parMultRC (Tpar *m);

#define perform(i1,j1,i2,j2,k)						\
   mult[k].treeA = treeA; mult[k].treeB = treeB;			\
   mult[k].level = level-1; mult[k].dlevel = dlevel+1; mult[k].plevel = plevel-1;\
   mult[k].res = &part[k];						\
   if (valA(i1) && valB(j1))						\
      { if (valA(i2) && valB(j2))					\
           { mult[k].nodeA1 = childA[i1]; mult[k].nodeB1 = childB[j1];	\
             mult[k].nodeA2 = childA[i2]; mult[k].nodeB2 = childB[j2];	\
             if (!PAR(plevel,dlevel) ||					\
                 pthread_create(&thr[k],NULL,(void*(*)(void*))parPerform,&mult[k]))	\
                { parPerform(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
        else 								\
           { mult[k].nodeA = childA[i1]; mult[k].nodeB = childB[j1];	\
             if (!PAR(plevel,dlevel) ||					\
	         pthread_create(&thr[k],NULL,(void*(*)(void*))parMult,&mult[k]))		\
                { parMult(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
      }									\
   else									\
      { if (valA(i2) && valB(j2))					\
           { mult[k].nodeA = childA[i2]; mult[k].nodeB = childB[j2];	\
             if (!PAR(plevel,dlevel) ||					\
	         pthread_create(&thr[k],NULL,(void*(*)(void*))parMult,&mult[k]))		\
                { parMult(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
	else { part[k] = empty; thr[k] = 0; }				\
      }									\

#define performRC(i1,j1,i2,j2,k)					\
   mult[k].treeA = treeA; mult[k].treeB = treeB;			\
   mult[k].level = level-1; mult[k].dlevel = dlevel+1; mult[k].plevel = plevel-1;\
   mult[k].row = row; mult[k].col = col;				\
   mult[k].res = &part[k];						\
   if (valA(i1) && valB(j1))						\
      { if (valA(i2) && valB(j2))					\
           { mult[k].nodeA1 = childA[i1]; mult[k].nodeB1 = childB[j1];	\
             mult[k].nodeA2 = childA[i2]; mult[k].nodeB2 = childB[j2];	\
             if (!PAR(plevel,dlevel) ||					\
		 pthread_create(&thr[k],NULL,(void*(*)(void*))parPerformRC,&mult[k]))	\
                { parPerformRC(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
        else 								\
           { mult[k].nodeA = childA[i1]; mult[k].nodeB = childB[j1];	\
             if (!PAR(plevel,dlevel) ||					\
	         pthread_create(&thr[k],NULL,(void*(*)(void*))parMultRC,&mult[k]))		\
                { parMultRC(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
      }									\
   else									\
      { if (valA(i2) && valB(j2))					\
           { mult[k].nodeA = childA[i2]; mult[k].nodeB = childB[j2];	\
             if (!PAR(plevel,dlevel) ||					\
	         pthread_create(&thr[k],NULL,(void*(*)(void*))parMultRC,&mult[k]))		\
                { parMultRC(&mult[k]); thr[k] = 0;			\
	        } 							\
	   }								\
	else { part[k] = empty;	thr[k] = 0; }				\
      }									\

static void parMult (Tpar *m);
static void parMultRC (Tpar *m);
static void parPerform (Tpar *m);
static void parPerformRC (Tpar *m);

static partition k2multRC (k2tree treeA, k2node nodeA, k2tree treeB, 
		k2node nodeB, uint level, uint64_t row, uint64_t col,
		uint dlevel, uint plevel)

   {
       //fprintf(stdout, "Time %llu\n", (t2-time_t1));
#if TIMEOUT
       if(time_diff() > TIMEOUT) {
           //fprintf(stdout, "Time %llu\n", user_diff());
           return empty;
       }
#endif

     partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB,sig;
     uint v;
     uint64_t lim;
     Tpar mult[4];
     pthread_t thr[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

     lim = ((uint64_t)1)<<(level-1);
     if (row != fullSide)
        { if (row >= lim) { sigA &= ~0x3; row -= lim; }
          else sigA &= ~0xC;
        }
     if (col != fullSide)
        { if (col >= lim) { sigB &= ~0x5; col -= lim; }
          else sigB &= ~0xA;
        }

        // base case, level = 1
     if (level == 1)
        { sig =  ((valA(0) & valB(0)) | (valA(1) & valB(2))) +
                (((valA(0) & valB(1)) | (valA(1) & valB(3))) << 1) +
                (((valA(2) & valB(0)) | (valA(3) & valB(2))) << 2) +
                (((valA(2) & valB(1)) | (valA(3) & valB(3))) << 3);
          return single (sig);
        }

	// combine them to form the 4 quadrants (now run in parallel)
     performRC(0,0,1,2,0);
     performRC(0,1,1,3,1);
     performRC(2,0,3,2,2);
     performRC(2,1,3,3,3);
     if (PAR(plevel,dlevel))
        for (v=0;v<4;v++)
	    if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static partition k2mult (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		         uint level, uint dlevel, uint plevel)

   {

       //fprintf(stdout, "Time %f\n", time_diff());
#if TIMEOUT
       if(time_diff() > TIMEOUT) {
           //fprintf(stdout, "Exit %f\n", time_diff());
           return empty;
       }
#endif

     partition part[4];
     Tpar mult[4];
     k2node childA[4],childB[4];
     uint sigA,sigB,sig;
     uint v;
     pthread_t thr[4];

     sigA = k2fillMappedChildren(treeA,nodeA,mapA,childA);
     sigB = k2fillMappedChildren(treeB,nodeB,mapB,childB);

        // base case, level = 1
     if (level == 1)
        { sig =  ((valA(0) & valB(0)) | (valA(1) & valB(2))) +
                (((valA(0) & valB(1)) | (valA(1) & valB(3))) << 1) +
                (((valA(2) & valB(0)) | (valA(3) & valB(2))) << 2) +
                (((valA(2) & valB(1)) | (valA(3) & valB(3))) << 3);
          return single (sig);
        }

	// compute the products of all the 4x4 combinations, now in parallel
     perform(0,0,1,2,0);
     perform(0,1,1,3,1);
     perform(2,0,3,2,2);
     perform(2,1,3,3,3);
     if (PAR(plevel,dlevel))
        for (v=0;v<4;v++)
	    if (thr[v]) pthread_join(thr[v],NULL);
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

static void parMult (Tpar *m)
   { *(m->res) = k2mult(m->treeA,m->nodeA,m->treeB,m->nodeB,
			m->level,m->dlevel,m->plevel);
   }

static void parMultRC (Tpar *m)
   { *(m->res) = k2multRC (m->treeA,m->nodeA,m->treeB,m->nodeB,
			   m->level,m->row,m->col,m->dlevel,m->plevel);
   }

static void parPerform (Tpar *m)

   { pthread_t thr;
     partition part1,part2;
     m->nodeA = m->nodeA1; m->nodeB = m->nodeB1;
     if (!PAR(m->plevel,m->dlevel) || pthread_create(&thr,NULL,(void*(*)(void*))parMult,m))
	{ parMult(m); thr = 0;
	}
     part2 = k2mult (m->treeA,m->nodeA2,m->treeB,m->nodeB2,m->level,m->dlevel,m->plevel);
     if (thr) pthread_join(thr,NULL);
     if (part2.tree != NULL)
        { if (m->res->tree == NULL) *(m->res) = part2;
	  else { part1 = *(m->res);
		 m->res->levels = myalloc(m->level*sizeof(uint64_t));
                 m->res->tree = k2merge(part1.tree,part1.len,
				    part2.tree,part2.len,m->level,
		                    &m->res->len,&m->res->elems,m->res->levels);
                 myfree(part1.tree); myfree(part1.levels);
                 myfree(part2.tree); myfree(part2.levels);
	       }
	}
   }

static void parPerformRC (Tpar *m)

   { pthread_t thr;
     partition part1,part2;
     m->nodeA = m->nodeA1; m->nodeB = m->nodeB1;
     if (!PAR(m->plevel,m->dlevel) || pthread_create(&thr,NULL,(void*(*)(void*))parMultRC,m))
	{ parMultRC(m); thr = 0;
	}
     part2 = k2multRC (m->treeA,m->nodeA2,m->treeB,m->nodeB2,
		       m->level,m->row,m->col,m->dlevel,m->plevel);
     if (thr) pthread_join(thr,NULL);
    if (part2.tree != NULL)
    { if (m->res->tree == NULL) *(m->res) = part2;
    else {
        part1 = *(m->res);
        m->res->levels = myalloc(m->level*sizeof(uint64_t));
        m->res->tree = k2merge(part1.tree,part1.len,part2.tree,part2.len,m->level,
                               &m->res->len,&m->res->elems,m->res->levels);
        myfree(part1.tree); myfree(part1.levels);
        myfree(part2.tree); myfree(part2.levels);
    }
    }
   }

        // (boolean) product of two matrices, assumed to be of the same side
        // only row of A and col of B are considered if not fullSide

matrix matMult (matrix A, matrix B)

   { return matMult1 (fullSide,A,B,fullSide);
   }

matrix matMult1 (uint64_t row, matrix A, matrix B, uint64_t col)

   { partition mult;
     matrix M;
     uint64_t hA,wB;

     DLEV = 2; OLEV = 1;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: product of matrices of different side\n");
          exit(1);
        }
     M = (matrix)myalloc(sizeof(struct s_matrix));
     M->logside = A->logside;
     matDims(A,NULL,NULL,&hA); matDims(B,NULL,&wB,NULL);
     M->width = wB; M->height = hA; 
     M->transposed = 0;
     if ((A->elems == 0) || (B->elems == 0))
	{ M->elems = 0; M->tree = NULL; }
     else
        { mapA = A->transposed ? mapTr : mapId;
          mapB = B->transposed ? mapTr : mapId;
	  if ((row == fullSide) && (col == fullSide))
	     mult = k2mult (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		            k2levels(A->tree),0,k2levels(A->tree));
	  else 
	     mult = k2multRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		              k2levels(A->tree),row,col,0,k2levels(A->tree));
          myfree (mult.levels);
          M->elems = mult.elems;
          if (M->elems == 0) M->tree = NULL;
          else M->tree = k2createFrom (k2levels(A->tree),mult.len,mult.tree,1);
	}
     return M;
   }

        // transitive closure of a matrix, pos says if it's + rather than *

	// k2sumRC with no col/row restrictions nor transposition, for speed 

static partition k2sum (k2tree treeA,k2node nodeA,k2tree treeB,k2node nodeB,
			uint level, uint has)

   { partition part[4];
     k2node childA[4],childB[4];
     uint sigA,sigB;
     uint v;

     if (has & 1) sigA = k2fillChildren(treeA,nodeA,childA);
     else sigA = 0;
     if (has & 2) sigB = k2fillChildren(treeB,nodeB,childB);
     else sigB = 0;

        // base case, level = 1
     if (level == 1) return single (sigA | sigB);

     for (v=0;v<4;v++) 
	 { has = valA(v) + 2*valB(v);
           if (has) 
              part[v] = k2sum (treeA,childA[v],treeB,childB[v],level-1,has);
	   else part[v] = empty;
  	 }
	// combine the 4 quadrants into a single matrix
     return compose(part,level,1);
   }

	// positive closure, add Id later to have reflexive

static void parMultMerge (Tpar *m)

   { *(m->res) = k2mult(m->treeA,m->nodeA,m->treeB,m->nodeB,m->level,m->dlevel,m->plevel);
     if (m->res->elems == 0) *(m->tres) = NULL; // ie AB
     else // add AB
	{
         k2tree taux = k2createFrom(m->level-1,m->res->len,m->res->tree,1);
	  *(m->res) = k2sum(m->treeB,m->nodeB,taux,k2root(taux),m->level,3);
	  k2destroy(taux);
	  *(m->tres) = k2createFrom(m->level,m->res->len,m->res->tree,1);
	}
   }

static void parMultMergeP (Tpar *m)

   { *(m->res) = k2mult(m->treeA,m->nodeA,m->treeB,m->nodeB,m->level,m->dlevel,m->plevel);
     if (m->res->elems == 0) *(m->tres) = NULL; // ie ABp
     else // add ABp
	{ k2tree taux = k2createFrom(m->level-1,m->res->len,m->res->tree,1);
	  *(m->res) = k2sum(m->treeA,m->nodeA,taux,k2root(taux),m->level,3);
	  k2destroy(taux);
	  *(m->tres) = k2createFrom(m->level,m->res->len,m->res->tree,1);
	}
   }

static partition k2clos (k2tree tree, k2node node, uint level, uint dlevel,
			 uint plevel)

   { k2tree tree1[4]; // AAp, ABp, BAp, BBp
     partition part1[4];
     k2tree tree2[4]; // AAs, ABs, BAs, BBs
     partition part2[4];
     partition part[4];
     partition answer;
     k2node child[4];
     uint v,sig,sigs;
     k2tree taux;
     uint share;
     pthread_t thr;

     sig = k2fillChildren(tree,node,child);

        // base case, level = 1
     if (level == 1)
        { sigs = (val(0) | (val(1) & val(2))) +
                 (val(1) << 1) +
                 (val(2) << 2) +
                 ((val(3) | (val(2) & val(1))) << 3);
          return single (sigs);
        }

     share = 0; // bit v of share = tree1[v] is shared by tree2[v]

	// first stage: nodes of AA can be intermediate
	// AAp = AA+
     if (!val(0)) { part1[0] = empty; tree1[0] = NULL; }
     else { part1[0] = k2clos (tree,child[0],level-1,dlevel+1,plevel);
	    tree1[0] = k2createFrom(level-1,part1[0].len,part1[0].tree,1);
	  }
	// ABp = AB | AAp x AB
	// tree1[1] == NULL means that ABp = tree1[1] = child 1 of tree = AB
     thr = 0;
     if ((tree1[0] == NULL) || !val(1))  // AAp or AB are empty
	  { part1[1] = empty; tree1[1] = NULL; } // ie part1[1] = AB
     else { Tpar mult;
	    mult.treeA = tree1[0]; mult.treeB = tree;
	    mult.nodeA = k2root(tree1[0]); mult.nodeB = child[1];
	    mult.level = level-1; mult.dlevel = dlevel+1; mult.plevel = plevel-1;
	    mult.res = &part1[1]; mult.tres = &tree1[1];
	    if (!PAR(plevel,dlevel) ||
		pthread_create(&thr,NULL,(void*(*)(void*))parMultMerge,&mult))
	       { parMultMerge(&mult); thr = 0; }
	  }
	// BAp = BA | BA x AAp
	// tree1[2] == NULL means that BAp = tree1[2] = child 2 of tree = BA
     if (!val(2) || (tree1[0] == NULL)) // ie BA or AAp are empty
	  { part1[2] = empty; tree1[2] = NULL; } // ie part1[2] = BA
     else { part1[2]= k2mult(tree,child[2],tree1[0],k2root(tree1[0]),
		     	     level-1,dlevel,plevel-1);
	    if (part1[2].elems == 0) tree1[2] = NULL; // ie BA
	    else // add BA
	       { if (val(2))
	            { taux = k2createFrom(level-1,part1[2].len,part1[2].tree,1);
		      part1[2]=k2sum(tree,child[2],taux,k2root(taux),level-1,3);
		      k2destroy(taux);
	            }
	         tree1[2] = k2createFrom(level-1,part1[2].len,part1[2].tree,1);
	       }
	  }
	// BBp = BB | BAp x AB (or also BB | BA x ABp)
	// tree1[3] == NULL means that BBp = tree1[3] = child 3 of tree = BB
     if (((tree1[2] == NULL) && !val(2)) || !val(1)) // BAp or AB are empty
	  { part1[3] = empty; tree1[3] = NULL; } // then BBp = BB
     else { if (tree1[2] == NULL) // means tree1[2] = BA
	       part1[3] = k2mult(tree,child[2],tree,child[1],level-1,dlevel,plevel);
	    else
	       part1[3] = k2mult(tree1[2],k2root(tree1[2]),tree,child[1],
			         level-1,dlevel,plevel);
	    if (part1[3].elems == 0) tree1[3] = NULL; // ie BB
	    else // add BB
	       { if (val(3))
		    { taux = k2createFrom(level-1,part1[3].len,part1[3].tree,1);
		      part1[3]=k2sum(tree,child[3],taux,k2root(taux),level-1,3);
		      k2destroy(taux);
		    }
		 tree1[3] = k2createFrom(level-1,part1[3].len,part1[3].tree,1);
	       }
	  }
	// closing computation of ABp
     if (thr) pthread_join(thr,NULL);

	// second stage: nodes of BB can also be intermediate
	// BBs = BBp+
     if (tree1[3]) part2[3]=k2clos (tree1[3],k2root(tree1[3]),level-1,dlevel+1,plevel);
     else if (val(3)) part2[3] = k2clos (tree,child[3],level-1,dlevel+1,plevel);
     else part2[3] = empty;
     if (part2[3].elems == 0) tree2[3] = NULL;
     else tree2[3] = k2createFrom(level-1,part2[3].len,part2[3].tree,1);
	// ABs = ABp | ABp x BBs
	// tree2[1] == NULL means that ABs = tree2[1] = ABp
	// but we then materialize ABs in this case too
     thr = 0;
     if (((tree1[1]==NULL) && !val(1)) || (tree2[3]==NULL)) // ABp or BBs empty
	  { part2[1] = empty; tree2[1] = NULL;  // then ABs = ABp
          }
     else { Tpar mult;
	    if (tree1[1] == NULL) // ie ABp = AB
	         { mult.treeA = tree; mult.nodeA = child[1]; }
	    else { mult.treeA = tree1[1]; mult.nodeA = k2root(tree1[1]); }
	    mult.treeB = tree2[3]; mult.nodeB = k2root(tree2[3]);
	    mult.level = level-1; mult.dlevel = dlevel+1; mult.plevel = plevel-1;
	    mult.res = &part2[1]; mult.tres = &tree2[1];
	    if (!PAR(plevel,dlevel) ||
		pthread_create(&thr,NULL,(void*(*)(void*))parMultMergeP,&mult))
	       { parMultMergeP(&mult); thr = 0; }
	  }
	// BAs = BAp | BBs x BAp
	// tree2[2] == NULL means that BAs = tree2[2] = BAp
	// but we then materialize BAs in this case too
     if ((tree2[3]==NULL) || ((tree1[2]==NULL) && !val(2))) // BAp or BBs empty
	  { part2[2] = empty; tree2[2] = NULL;  // then BAs = BAp
          }
     else { if (tree1[2] == NULL) // ie BAp = BA
	         part2[2] = k2mult(tree2[3],k2root(tree2[3]),
			 	   tree,child[2],level-1,dlevel,plevel-1);
	    else part2[2] = k2mult(tree2[3],k2root(tree2[3]),
			           tree1[2],k2root(tree1[2]),level-1,dlevel,plevel-1);
	    if (part2[2].elems == 0) tree2[2] = NULL; // ie BAp
	    else // add BAp
	       { taux = k2createFrom(level-1,part2[2].len,part2[2].tree,1);
	         if (tree1[2] == NULL) // ie BAp = BA
		    part2[2]=k2sum(tree,child[2],taux,k2root(taux),level-1,3);
		 else
		    part2[2]=k2sum(tree1[2],k2root(tree1[2]),taux,k2root(taux),
				   level-1,3);
		 k2destroy(taux);
	         tree2[2] = k2createFrom(level-1,part2[2].len,part2[2].tree,1);
	       }
          }
     if (part2[2].elems == 0) // tree2[2]==NULL also holds, materialize BAs
        { if (tree1[2] == NULL) // BAp = BA
	     { if (val(2))
		  { part2[2] = k2extract(tree,child[2],level-1);
	            tree2[2]=k2createFrom(level-1,part2[2].len,part2[2].tree,1);
		  }
	        // else actually empty
	     }
	       // avoid copying:
	       // part2[2] = k2extract(tree1[2],k2root(tree[2]),level-1);
          else
	     { share |= 0x4;  // mark that we are sharing the tree
	       part2[2] = part1[2]; tree2[2] = tree1[2];
	     }
	}
	// closing computation of ABs
     if (thr) pthread_join(thr,NULL);
     if (part2[1].elems == 0) // tree2[1]==NULL also holds, materialize ABs
        { if (tree1[1] == NULL) // ABp = AB
	     { if (val(1))
		  { part2[1] = k2extract(tree,child[1],level-1);
	            tree2[1]=k2createFrom(level-1,part2[1].len,part2[1].tree,1);
		  }
	        // else actually empty
	     }
	       // avoid a copy:
	       // part2[1] = k2extract(tree1[1],k2root(tree[1]),level-1);
          else
	     { share |= 0x2;  // mark that we are sharing the tree
	       part2[1] = part1[1]; tree2[1] = tree1[1];
	     }
	}
	// AAs = AAp | ABp x BAs (or also AAp | ABs x BAp)
	// tree2[0] == NULL means that AAs = tree2[0] = AAp
	// but we then materialize AAs in this case too
     if (((tree1[1]==NULL) && !val(1)) || (tree2[2]==NULL)) // ABp or BAs empty
	{ part2[0] = empty; tree2[0] = NULL;  // then AAs = AAp
        }
     else
	{ if (tree1[1] == NULL) // ABp = AB
	     part2[0] = k2mult(tree,child[1],tree2[2],k2root(tree2[2]),
			       level-1,dlevel,plevel);
	  else
	     part2[0] = k2mult(tree1[1],k2root(tree1[1]),
			       tree2[2],k2root(tree2[2]),level-1,dlevel,plevel);
	  if (part2[0].elems == 0) tree2[0] = NULL; // ie AAp
          else // add AAp
	     { if (tree1[0] != NULL)
	          { taux = k2createFrom(level-1,part2[0].len,part2[0].tree,1);
	            part2[0] = k2sum(tree1[0],k2root(tree1[0]),
				     taux,k2root(taux),level-1,3);
	            k2destroy(taux);
	          }
	       tree2[0] = k2createFrom(level-1,part2[0].len,part2[0].tree,1);
	     }
	}
     if (tree2[0] == NULL) // then AAs = AAp
	{ share |= 0x1; // this time could just set tree1[0] & part1[0] empty
	  part2[0] = part1[0]; tree2[0] = tree1[0];
	}

	// combine the 4 quadrants into a single matrix
     answer = compose(part2,level,0);

     for (v=0;v<4;v++)
	 { if (!(share & (1<<v)))
	      { if (tree1[v] != NULL) k2destroy(tree1[v]);
	        myfree (part1[v].levels);
	      }
	   if (tree2[v] != NULL) k2destroy(tree2[v]);
	   myfree (part2[v].levels);
	 }

     return answer;
   }

matrix matClos (matrix A, uint pos)

   { matrix M,S,Id;
     partition aux;
     uint levels;
     uint side = mmax(A->width,A->height);

     //DLEV = 6; OLEV = 5;
     DLEV = 2; OLEV = 1;
     if (A->elems == 0)
        if (pos) return matEmpty(side,side);
	else return matId(side);

     mapA = mapB = mapId; // just for k2mult and k2extract, not reentrant!

     levels = k2levels(A->tree);
     aux = k2clos(A->tree,k2root(A->tree),levels,0,levels);
     myfree(aux.levels);
     M = (matrix)myalloc(sizeof(struct s_matrix));
     M->width = M->height = side;
     M->logside = A->logside; M->transposed = A->transposed;
     M->elems = aux.elems;
     M->tree = k2createFrom(levels,aux.len,aux.tree,1);
	// add Id if not pos
     if (!pos)
        { Id = matId(M->width);
          S = M; M = matSum(M,Id); matDestroy(S);
          matDestroy(Id);
	}
     return M;
   }

matrix matClos0 (matrix A, uint pos)

   { matrix M,P,S,Id;
     uint64_t elems;

     if (!pos)
        { Id = matId(mmax(A->width,A->height));
          A = matSum(A,Id);
          matDestroy(Id);
	}
     elems = A->elems;
     if (elems == 0) return matCopy(A); // can only be pos, A not to destroy
     P = matMult (A,A);
     S = matSum (A,P);
     if (!pos) matDestroy(A);
     while (S->elems != elems)
	{ elems = S->elems;
	  matDestroy(P);
	  P = matMult(S,S);
	  M = matSum(S,P);
	  matDestroy(S);
	  S = M;
	}
     matDestroy(P);
     return S;
   }

        // versions to chose row row or column col at the end
	// here we restrict first to the row/column even if we lose the
	// exponential increase in the paths

	// if coltest not null, checks that the cell *coltest exists and 
	// returns there 1 or 0 as soon as it can
static matrix matClosRow (uint64_t row, matrix I, matrix A, uint pos, 
			  uint64_t *coltest)

   { matrix M,P,S,E;
     uint64_t elems;
     uint64_t dim = mmax(A->height,A->width);

     if (I == NULL)
        { if (pos) E = matEmpty(dim,dim);
          else E = matOne(dim,dim,row,row);
          S = matOr1 (row,A,E,fullSide);
          matDestroy(E);
        }
     else
        { if (pos) S = matMult1(row,I,A,fullSide);
          else { E = matEmpty(dim,dim);
                 S = matOr1(row,I,E,fullSide);
                 matDestroy(E);
               }
        }
     elems = S->elems;
     if (coltest && matAccess(S,row,*coltest))
	{ matDestroy(S); *coltest = 1; return NULL; }
     P = matMult (S,A);
     M = S; S = matOr (S,P); matDestroy(M);
     while (S->elems != elems)
	{ if (coltest && matAccess(S,row,*coltest))
	     { matDestroy(S); matDestroy(P); *coltest = 1; return NULL; }
	  elems = S->elems;
	  M = P; P = matMult(P,A); matDestroy(M);
	  M = S; S = matOr(S,P); matDestroy(M);
	}
     matDestroy(P);
     if (coltest) { matDestroy(S); *coltest = 0; return NULL; }
     return S;
   }

        // versions to choose one row or one column, or both 

matrix matClos1 (uint64_t row, matrix A, uint pos, uint64_t col)

   { uint64_t nrow,ncol;
     uint64_t side = mmax(A->width,A->height);
     uint64_t cell[2];
     uint64_t test;
     struct s_matrix At; 
     matrix M;
     if (row == fullSide)
	{ if (col == fullSide) return matClos(A,pos);
	  else 
             { At = matTranspose(A);
               M = matClosRow(col,NULL,&At,pos,NULL);
	       *M = matTranspose(M);
               return M;
             }
	}
     else
        { if (col == fullSide) return matClosRow(row,NULL,A,pos,NULL);
	}
	// both row and col
     nrow = matCollect (A,0,fullSide,col,col,NULL);
     ncol = matCollect (A,row,row,0,fullSide,NULL);
     if (ncol < nrow)
	{ test = row;
	  At = matTranspose(A);
	  matClosRow(col,NULL,&At,pos,&test); // does not return matrix
	}
     else
	{ test = col;
	  matClosRow(row,NULL,A,pos,&test); // does not return matrix
	}
     if (test) return matOne (side, side, row, col);
     else return matEmpty(side,side);
   }

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)

matrix matMultClos1 (uint64_t row, matrix A, matrix B, uint pos, uint64_t col)

   { uint64_t nrow,ncol;
     uint64_t side = mmax(A->width,A->height);
     uint64_t test;
     matrix M,M1;
     if (row == fullSide)
        { if (col == fullSide)  // nothing special
             { M1 = matClos(B,pos);
               M = matMult(A,M1);
               matDestroy(M1);
               return M;
             }
          else // it's much better to start with restricted B
             { M1 = matClos1(fullSide,B,pos,col);
               M = matMult(A,M1);
               matDestroy(M1);
               return M;
             }
        }
     else
        { if (col == fullSide) return matClosRow(row,A,B,pos,NULL);
        }
        // both row and col
     test = col;
     matClosRow(row,A,B,pos,&test); // returns no matrix, just test
     if (test) return matOne (side,side,row,col);
     return matEmpty(side,side);
   }

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)

matrix matClosMult1 (uint64_t row, matrix A, uint pos, matrix B, uint64_t col)

   { matrix M;
     struct s_matrix At,Bt;
     At = matTranspose(A);
     Bt = matTranspose(B);
     M = matMultClos1(col,&Bt,&At,pos,row);
     *M = matTranspose(M);
     return M;
   }

        // multiplies M by vector V. allocates and returns the result M x V
        // use matrix transposition to do V^T x M = M^T x V (transposed vector)
        // V is assumed to be of size A->width, the output is of size A->height

typedef struct
   { k2tree tree;
     k2node node;
     uint level,dlevel;
     double *V,*W;
   } Tparv;

static void parVecMult (Tparv *m);

static void k2vectorMult (k2tree tree, k2node node, uint level,
		          double *V, double *W, uint dlevel)

   { uint v;
     uint lshift;
     uint sig;
     k2node child[4];
     Tparv m;
     pthread_t thr;

     if (level == 0)
	{ *W += *V; return; }
     level--; dlevel++; lshift = ((uint)1) << level;
     sig = k2fillMappedChildren(tree,node,mapA,child);
     if (!PAR(level,dlevel))
        { for (v=0;v<4;v++)
              if (val(v))
                 k2vectorMult(tree,child[v],level,
                          V + lshift * (v & 1), W + lshift * (v >> 1),dlevel);
	  return;
	}
	// parallelize, but not the updates to W
	// 0 and 2 can go in parallel, but only if both have to run
     thr = 0;
     m.tree = tree; m.level = level; m.dlevel = dlevel;
     m.V = V; m.W = W;
     if (val(0))
	{ m.node = child[0];
          if (!val(2) || pthread_create(&thr,NULL,(void*(*)(void*))parVecMult,&m))
	     { parVecMult(&m); thr = 0; }
        }
     if (val(2))
	k2vectorMult(tree,child[2],level,V,W+lshift,dlevel);
     if (thr) pthread_join(thr,NULL);
	// 1 and 3 can go in parallel, but only if both have to run
     thr = 0;
     if (val(1))
	{ m.node = child[1];
	  m.V += lshift;
	  if (!val(3) || pthread_create(&thr,NULL,(void*(*)(void*))parVecMult,&m))
	     { parVecMult(&m); thr = 0; }
	}
     if (val(3))
	k2vectorMult(tree,child[3],level,V+lshift,W+lshift,dlevel);
     if (thr) pthread_join(thr,NULL);
   }

static void parVecMult (Tparv *m)

   { k2vectorMult(m->tree,m->node,m->level,m->V,m->W,m->dlevel);
   }

double *matVectorMult (matrix A, double *V)

   { double *W; // W = MxV
     uint i;

     DLEV = 3; OLEV = 10;  // does not outperform sequential :-(

     if (A->height == 0) return NULL;
     W = (double*)myalloc(A->height * sizeof(double));
     for (i=0;i<A->height;i++) W[i] = 0;
     if (A->elems)
        { mapA = A->transposed ? mapTr : mapId;
          k2vectorMult(A->tree,k2root(A->tree),k2levels(A->tree),V,W,0);
	}
     return W;
   }

