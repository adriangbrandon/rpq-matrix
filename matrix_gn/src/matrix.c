
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"

typedef struct {
   k2node nodeA,nodeB;
   uint dist;
   } queue;

static queue *Q = NULL;
static uint64_t head,tail,size;

	// creates matrix of width x height with n cells (2n ints row,col) 
	// reorders cells array

static uint numbits (uint64_t n)

   { uint bits = 0;
     while (n)
	{ n = n>>1; bits++; }
     return bits ? bits : 1;
   }
	
matrix matCreate (uint64_t height, uint64_t width, uint64_t n, uint64_t *cells)

   { matrix mat = (matrix)malloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(max(width,height)-1);
     mat->elems = n;
     mat->transposed = 0;
     if (n == 0) mat->tree = NULL;
     else mat->tree = k2create(n,mat->logside,cells);
     return mat;
   }
     
matrix matCreate32 (uint64_t height, uint64_t width, uint64_t n, uint32_t *cells)

   { matrix mat = (matrix)malloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(max(width,height)-1);
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
     free(M);
   }

        // creates an empty matrix
matrix matEmpty (uint64_t height, uint64_t width)

   { matrix mat = (matrix)malloc(sizeof(struct s_matrix));
     mat->width = width ? width : 1;
     mat->height = height ? height : 1;
     mat->logside = numbits(max(width,height)-1);
     mat->elems = 0;
     mat->transposed = 0;
     mat->tree = NULL;
     return mat;
   }

        // creates an identity matrix

matrix matId (uint64_t side)

   { matrix mat = (matrix)malloc(sizeof(struct s_matrix));
     uint64_t *data;
     uint64_t nodes;
     side = side ? side : 1;
     mat->width = side;
     mat->height = side;
     mat->logside = numbits(side-1);
     mat->elems = side;
     mat->transposed = 0;
     nodes = (((uint64_t)1) << mat->logside)-1;
     data = (uint64_t*)malloc(((4*nodes+w-1)/w)*sizeof(uint64_t));
     memset(data,0x99,(nodes+1)/2);
     mat->tree = k2createFrom(mat->logside,4*nodes,data,1);
     return mat;
   }
     
	// creates a new copy of A, with its own data

matrix matCopy (matrix A)

   { matrix M = (matrix)malloc(sizeof(struct s_matrix));
     *M = *A;
     if (A->tree != NULL) M->tree = k2copy(A->tree);
     return M;
   }

	// transposes matrix M

void matTranspose (matrix M)

   { M->transposed = ~M->transposed;
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
      M = (matrix)malloc(sizeof(struct s_matrix));
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
     if (r2 == fullSide) r2 = M->height-1; // already in concrete repres
     if (c2 == fullSide) c2 = M->width-1;
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
     else { M = (matrix)malloc(sizeof(struct s_matrix));
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
	{ M->width = max(A->width,B->width);
     	  M->height = max(A->height,B->height);
	}
     else if (A->elems > B->elems) // the transposition of A dominates
	{ M->width = max(A->width,B->height);
     	  M->height = max(A->height,B->width);
	}
     else // A->elems <= B->elems, the transposition of B dominates
	{ M->width = max(B->width,A->height);
     	  M->height = max(B->height,A->width);
	}
     return M;
   }

        // version with one row or one column, or both

typedef struct {
   short has;
   uint64_t lim;
   k2node nodeA;
   k2node nodeB;
   uint64_t row,col;
   } queue4;

static uint mapId[] = { 0, 1, 2, 3 };
static uint mapTr[] = { 0, 2, 1, 3 };

static uint *mapA,*mapB;

static uint64_t *k2sum1 (k2tree A, k2tree B, uint64_t row, uint64_t col,
                         uint64_t *len, uint64_t *telems)

   { uint64_t *treeM;
     uint64_t ptr;
     uint level,v;
     uint64_t elems,size,lenA,lenB,lim;
     queue4 *Q;

     level = k2levels(B);
     lenA = bitsLength(k2bits(B));
     lenB = bitsLength(k2bits(B));

     treeM = (uint64_t*)malloc(((lenA+lenB+w-1)/w)*sizeof(uint64_t));
     size = lenA + lenB + 1;
     Q = (queue4*)malloc(size*sizeof(queue4));
     ptr = 0;
     elems = 0;
     head = 0; tail = 1;
     Q[0].has = (A != NULL) + 2 * (B != NULL);
     Q[0].lim = ((uint64_t)1)<<level;
     Q[0].nodeA = k2root(A);
     Q[0].nodeB = k2root(B);
     Q[0].row = Q[0].col = 0; 
     while (head != tail)
        { uint has,has1,has2;
          k2node nodeA,nodeB;
	  uint64_t lim = Q[head].lim >> 1;
	  uint64_t row0 = Q[head].row;
	  uint64_t col0 = Q[head].col;
	  uint64_t row1,col1;
          has = Q[head].has;
          nodeA = Q[head].nodeA;
          nodeB = Q[head].nodeB;
	  if ((has & 1) == 0) // then node of B not null
	     { for (v=0;v<4;v++)
	  	   { row1 = row0 + (v >> 1) * lim;
	  	     col1 = col0 + (v & 1) * lim;
	             has2 = ((row == fullSide) || ((row>=row1)&&(row<row1+lim)))
	                && ((col == fullSide) || ((col>=col1)&&(col<col1+lim)));
	             if (has2) has2 = k2hasChild(B,nodeB,mapB[v]);
	             bitsWriteA(treeM,ptr++,has2);
	             if (has2)
	                { if (lim > 1)
			     { Q[tail].lim = lim;
		               Q[tail].has = 2;
			       Q[tail].nodeB = k2child(B,nodeB,mapB[v]);
		               Q[tail].row = row1;
		     	       Q[tail].col = col1;
		               tail = (tail+1)%size;
		             }
			  else elems++;
		        }
		   }
	     }
	  else if ((has & 2) == 0) // node of A not null anyway
	     { for (v=0;v<4;v++)
	  	   { row1 = row0 + (v >> 1) * lim;
	  	     col1 = col0 + (v & 1) * lim;
	             has1 = ((row == fullSide) || ((row>=row1)&&(row<row1+lim)))
	                && ((col == fullSide) || ((col>=col1)&&(col<col1+lim)));
	             if (has1) has1 = k2hasChild(A,nodeA,mapA[v]);
	             bitsWriteA(treeM,ptr++,has1);
		     if (has1)
			{ if (lim > 1)
		             { Q[tail].has = 1;
			       Q[tail].lim = lim;
			       Q[tail].nodeA = k2child(A,nodeA,mapA[v]);
		               Q[tail].row = row1;
		     	       Q[tail].col = col1;
		               tail = (tail+1)%size;
		             }
			  else elems++;
			}
		   }
	     }
	  else // both A and B nodes exist
	     { for (v=0;v<4;v++)
	  	   { row1 = row0 + (v >> 1) * lim;
	  	     col1 = col0 + (v & 1) * lim;
	             has1 = has2 =
	                   ((row == fullSide) || ((row>=row1)&&(row<row1+lim)))
	                && ((col == fullSide) || ((col>=col1)&&(col<col1+lim)));
	             if (has1) has1 = k2hasChild(A,nodeA,mapA[v]);
	             if (has2) has2 = k2hasChild(B,nodeB,mapB[v]);
	             bitsWriteA(treeM,ptr++,has1|has2);
	             if (has1 | has2)
			{ if (lim > 1)
		             { Q[tail].has = 2*has2 + has1;
			       Q[tail].lim = lim;
			       if (has1) 
				  Q[tail].nodeA = k2child(A,nodeA,mapA[v]);
			       if (has2) 
				  Q[tail].nodeB = k2child(B,nodeB,mapB[v]);
		               Q[tail].row = row1;
		     	       Q[tail].col = col1;
		               tail = (tail+1)%size;
		             }
			  else elems++;
		        }
		   }
	     }
	  head = (head+1)%size;
        }
     free(Q); Q = NULL; // safety
     treeM = (uint64_t*)realloc(treeM,((ptr+w-1)/w)*sizeof(uint64_t));
     *telems = elems;
     *len = ptr;
     return treeM;
   }

matrix matSum1 (matrix A, matrix B, uint64_t row, uint64_t col)

   { uint64_t *sum;
     matrix M;
     uint64_t wA,hA,wB,hB;
     uint64_t len;

     if ((row == fullSide) && (col == fullSide)) return matSum(A,B);

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: sum of matrices of different side\n");
          exit(1);
        }
     M = (matrix)malloc(sizeof(struct s_matrix));
     M->logside = A->logside;
     matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
     M->width = max(wA,wB); M->height = max(hA,hB); 
     M->transposed = 0;
     mapA = A->transposed ? mapTr : mapId;
     mapB = B->transposed ? mapTr : mapId;
     sum = k2sum1 (A->tree,B->tree,row,col,&len,&M->elems);
     if (M->elems == 0) M->tree = NULL;
     else M->tree = k2createFrom (k2levels(A->tree),len,sum,1);
     return M;
   }

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

	// (boolean) product of two matrices, assumed to be of the same side

#define prepare(i1,j1,i2,j2)						\
     if (valA[i1] && valB[j1]) 						\
 	  part1 = k2mult (treeA,childA[i1],treeB,childB[j1],level-1);	\
     else part1 = empty;						\
     if (valA[i2] && valB[j2]) 						\
 	  part2 = k2mult (treeA,childA[i2],treeB,childB[j2],level-1);	\
     else part2 = empty;						\

#define prepareRC(i1,j1,i2,j2)						\
     if (valA[i1] && valB[j1]) 						\
 	  part1 = k2multRC (treeA,childA[i1],treeB,childB[j1],level-1,row,col);\
     else part1 = empty;						\
     if (valA[i2] && valB[j2]) 						\
 	  part2 = k2multRC (treeA,childA[i2],treeB,childB[j2],level-1,row,col);\
     else part2 = empty;						\

#define combine(i1,j1,i2,j2,k)						\
     if (part1.tree == NULL) part[k] = part2; 				\
     else if (part2.tree == NULL) part[k] = part1;			\
     else								\
        { part[k].levels = malloc((level-1)*sizeof(uint64_t));		\
	  part[k].tree = k2merge(part1.tree,part1.len,part2.tree,part2.len,\
			   level-1,&part[k].len,&part[k].elems,part[k].levels);\
	  free(part1.tree); free(part1.levels);				\
	  free(part2.tree); free(part2.levels);				\
	}

typedef struct {
   uint64_t *tree,*levels;
   uint64_t elems,len;
   } partition;

static partition empty = { NULL, NULL, 0, 0 };

static partition k2multRC (k2tree treeA, k2node nodeA, k2tree treeB, 
		k2node nodeB, uint level, uint64_t row, uint64_t col)

   { partition part1,part2,part[4];
     k2node childA[4],childB[4];
     uint valA[4],valB[4];
     uint64_t ptrs[4];
     uint u,v,val,l,old;
     uint64_t ptr,lim;
     partition answer;

     for (v=0;v<4;v++) valA[v] = valB[v] = 0;
     lim = ((uint64_t)1)<<(level-1);
     if ((row == fullSide) || (row < lim))
        { valA[0] = k2hasChild(treeA,nodeA,mapA[0]);
          valA[1] = k2hasChild(treeA,nodeA,mapA[1]);
	}
     if ((row == fullSide) || (row >= lim))
	{ valA[2] = k2hasChild(treeA,nodeA,mapA[2]);
          valA[3] = k2hasChild(treeA,nodeA,mapA[3]);
	  if (row != fullSide) row -= lim;
	}
     if ((col == fullSide) || (col < lim))
        { valB[0] = k2hasChild(treeB,nodeB,mapB[0]);
          valB[2] = k2hasChild(treeB,nodeB,mapB[2]);
	}
     if ((col == fullSide) || (col >= lim))
	{ valB[1] = k2hasChild(treeB,nodeB,mapB[1]);
          valB[3] = k2hasChild(treeB,nodeB,mapB[3]);
	  if (col != fullSide) col -= lim;
	}
	// base case, level = 1
     if (level == 1)
	{ val =  ((valA[0] & valB[0]) | (valA[1] & valB[2])) + 
                (((valA[0] & valB[1]) | (valA[1] & valB[3])) << 1) +
                (((valA[2] & valB[0]) | (valA[3] & valB[2])) << 2) +
                (((valA[2] & valB[1]) | (valA[3] & valB[3])) << 3);
 	  if (val == 0) return empty;
	  answer.elems = popcount(val);
          answer.len = 4;
          answer.tree = (uint64_t*)malloc(((4+w-1)/w)*sizeof(uint64_t));
	  answer.tree[0] = val;
	  answer.levels = (uint64_t*)malloc(1*sizeof(uint64_t));
	  answer.levels[0] = 1;
	  return answer;
	}

	// compute the products of all the 4x4 combinations
     for (v=0;v<4;v++) 
	 { if (valA[v]) childA[v] = k2child(treeA,nodeA,mapA[v]);
	   if (valB[v]) childB[v] = k2child(treeB,nodeB,mapB[v]);
	 }
	// combine them to form the 4 quadrants
     prepareRC(0,0,1,2);
       combine(0,0,1,2,0);
     prepareRC(0,1,1,3);
       combine(0,1,1,3,1);
     prepareRC(2,0,3,2);
       combine(2,0,3,2,2);
     prepareRC(2,1,3,3);
       combine(2,1,3,3,3);
	// combine the 4 quadrants into a single matrix
     answer.elems = part[0].elems+part[1].elems+part[2].elems+part[3].elems;
     if (answer.elems == 0) answer = empty;
     else
       { answer.len = 4 + part[0].len+part[1].len+part[2].len+part[3].len;
         answer.levels = (uint64_t*)malloc(level*sizeof(uint64_t));
         answer.levels[level-1] = 1;
         answer.tree = (uint64_t*)malloc((1+(answer.len+w-1)/w)
				         *sizeof(uint64_t)); //+1 for copyBits
         for (v=0;v<4;v++) 
	     { bitsWriteA(answer.tree, v, part[v].elems);
	       ptrs[v] = 0;
	     }
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
        }
	// free all the components 
     for (u=0;u<4;u++)
	 { if (part[u].tree != NULL) free(part[u].tree);
	   if (part[u].levels != NULL) free (part[u].levels);
	 }
     return answer;
   }

static partition k2mult (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
		         uint level)

   { partition part1,part2,part[4];
     k2node childA[4],childB[4];
     uint valA[4],valB[4];
     uint64_t ptrs[4];
     uint u,v,val,l,old;
     uint64_t ptr;
     partition answer;

     for (v=0;v<4;v++) 
         { valA[v] = k2hasChild(treeA,nodeA,mapA[v]);
	   valB[v] = k2hasChild(treeB,nodeB,mapB[v]);
         }
	// base case, level = 1
     if (level == 1)
	{ val =  ((valA[0] & valB[0]) | (valA[1] & valB[2])) + 
                (((valA[0] & valB[1]) | (valA[1] & valB[3])) << 1) +
                (((valA[2] & valB[0]) | (valA[3] & valB[2])) << 2) +
                (((valA[2] & valB[1]) | (valA[3] & valB[3])) << 3);
 	  if (val == 0) return empty;
	  answer.elems = popcount(val);
          answer.len = 4;
          answer.tree = (uint64_t*)malloc(((4+w-1)/w)*sizeof(uint64_t));
	  answer.tree[0] = val;
	  answer.levels = (uint64_t*)malloc(1*sizeof(uint64_t));
	  answer.levels[0] = 1;
	  return answer;
	}

	// compute the products of all the 4x4 combinations
     for (v=0;v<4;v++) 
	 { if (valA[v]) childA[v] = k2child(treeA,nodeA,mapA[v]);
	   if (valB[v]) childB[v] = k2child(treeB,nodeB,mapB[v]);
	 }
	// combine them to form the 4 quadrants
     prepare(0,0,1,2);
     combine(0,0,1,2,0);
     prepare(0,1,1,3);
     combine(0,1,1,3,1);
     prepare(2,0,3,2);
     combine(2,0,3,2,2);
     prepare(2,1,3,3);
     combine(2,1,3,3,3);
	// combine the 4 quadrants into a single matrix
     answer.elems = part[0].elems+part[1].elems+part[2].elems+part[3].elems;
     if (answer.elems == 0) answer = empty;
     else
       { answer.len = 4 + part[0].len+part[1].len+part[2].len+part[3].len;
         answer.levels = (uint64_t*)malloc(level*sizeof(uint64_t));
         answer.levels[level-1] = 1;
         answer.tree = (uint64_t*)malloc((1+(answer.len+w-1)/w)
				         *sizeof(uint64_t)); //+1 for copyBits
         for (v=0;v<4;v++) 
	     { bitsWriteA(answer.tree, v, part[v].elems);
	       ptrs[v] = 0;
	     }
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
        }
	// free all the components 
     for (u=0;u<4;u++)
	 { if (part[u].tree != NULL) free(part[u].tree);
	   if (part[u].levels != NULL) free (part[u].levels);
	 }
     return answer;
   }

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide

matrix matMult (matrix A, matrix B)

   { return matMult1 (A,fullSide,B,fullSide);
   }

matrix matMult1 (matrix A, uint64_t rowA, matrix B, uint64_t colB)

   { partition mult;
     matrix M;
     uint64_t hA,wB;

     if (A->logside != B->logside)
        { fprintf(stderr,"Error: product of matrices of different side\n");
          exit(1);
        }
     M = (matrix)malloc(sizeof(struct s_matrix));
     M->logside = A->logside;
     matDims(A,NULL,NULL,&hA); matDims(B,NULL,&wB,NULL);
     M->width = wB; M->height = hA; 
     M->transposed = 0;
     if ((A->elems == 0) || (B->elems == 0))
	{ M->elems = 0; M->tree = NULL; }
     else
        { mapA = A->transposed ? mapTr : mapId;
          mapB = B->transposed ? mapTr : mapId;
	  if ((rowA == fullSide) && (colB == fullSide))
	     mult = k2mult (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		         k2levels(A->tree));
	  else 
	     mult = k2multRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
		         k2levels(A->tree),rowA,colB);
          if (mult.levels != NULL) free (mult.levels);
          M->elems = mult.elems;
          if (M->elems == 0) M->tree = NULL;
          else M->tree = k2createFrom (k2levels(A->tree),mult.len,mult.tree,1);
	}
     return M;
   }

        // transitive closure of a matrix, pos says if it's + rather than *

matrix matClos (matrix A, uint pos)

   { matrix M,P,S,Id;
     uint64_t elems;
     uint transp = A->transposed;

     A->transposed = 0; // may be slightly more cache-friendly
     if (!pos)
        { Id = matId(max(A->width,A->height));
          M = matSum(A,Id);
          matDestroy(Id);
          A = M;
	}
     elems = A->elems;
     if (elems == 0) return matCopy(A);
     P = matMult (A,A);
     S = matSum (A,P);
     while (S->elems != elems)
	{ elems = S->elems;
	  matDestroy(P);
	  P = matMult(S,S);
	  M = matSum(S,P);
	  matDestroy(S);
	  S = M;
	}
     matDestroy(P);
     S->transposed = transp;
     return S;
   }

        // versions to chose row row or column col at the end
	// here we restrict first to the row/column even if we lose the
	// exponential increase in the paths

	// if coltest not null, checks that the cell *coltest exists and 
	// returns there 1 or 0 as soon as it can
static matrix matClosRow (matrix A, uint pos, uint64_t row, uint64_t *coltest)

   { matrix Id,M,P,S,Ar;
     uint64_t elems;

     Id = matId(max(A->height,A->width));
     Ar = matMult1 (A,row,Id,fullSide);
     if (!pos)
	{ M = matMult1(Id,row,Id,row); // just the cell (row,row)
	  matDestroy(Id);
	  Id = M;
	  M = matSum(Id,Ar);
	  matDestroy(Ar);
	  Ar = M;
	}
     matDestroy(Id);
     elems = Ar->elems;
     if (coltest && matAccess(Ar,row,*coltest))
	{ matDestroy(Ar); *coltest = 1; return NULL; }
     P = matMult (Ar,A);
     S = matSum (P,Ar); 
     matDestroy(Ar);
     while (S->elems != elems)
	{ if (coltest && matAccess(S,row,*coltest))
	     { matDestroy(S); matDestroy(P); *coltest = 1; return NULL; }
	  elems = S->elems;
	  M = matMult(P,A);
	  matDestroy(P);
	  P = M;
	  M = matSum(S,P);
	  matDestroy(S);
	  S = M;
	}
     matDestroy(P);
     if (coltest) { *coltest = 0; return NULL; }
     return S;
   }

static matrix matClosCol (matrix A, uint pos, uint64_t col, uint64_t *rowtest)

   { matrix Id,M,P,S,Ac;
     uint64_t elems;

     Id = matId(max(A->height,A->width));
     Ac = matMult1 (Id,fullSide,A,col);
     if (!pos)
	{ M = matMult1(Id,col,Id,col); // just the cell (col,col)
	  matDestroy(Id);
	  Id = M;
	  M = matSum(Id,Ac);
	  matDestroy(Ac);
	  Ac = M;
	}
     matDestroy(Id);
     elems = Ac->elems;
     if (rowtest && matAccess(Ac,*rowtest,col))
	{ matDestroy(Ac); *rowtest = 1; return NULL; }
     P = matMult (A,Ac);
     S = matSum (P,Ac); 
     matDestroy(Ac);
     while (S->elems != elems)
	{ if (rowtest && matAccess(S,*rowtest,col))
	     { matDestroy(S); matDestroy(P); *rowtest = 1; return NULL; }
	  elems = S->elems;
	  M = matMult(A,P);
	  matDestroy(P);
	  P = M;
	  M = matSum(S,P);
	  matDestroy(S);
	  S = M;
	}
     matDestroy(P);
     if (rowtest) { *rowtest = 0; return NULL; }
     return S;
   }

        // versions to choose one row or one column, or both 

matrix matClos1 (matrix A, uint pos, uint64_t row, uint64_t col)

   { uint64_t nrow,ncol;
     uint64_t side = max(A->width,A->height);
     uint64_t cell[2];
     uint64_t test;
     if (row == fullSide)
	{ if (col == fullSide) return matClos(A,pos);
	  else return matClosCol(A,pos,col,NULL);
	}
     else
        { if (col == fullSide) return matClosRow(A,pos,row,NULL);
	}
	// both row and col
     nrow = matCollect (A,0,fullSide,col,col,NULL);
     ncol = matCollect (A,row,row,0,fullSide,NULL);
     if (ncol < nrow)
	{ test = row;
	  matClosCol(A,pos,col,&test);
	}
     else
	{ test = col;
	  matClosRow(A,pos,row,&test);
	}
     if (test)
	{ cell[0] = row; cell[1] = col;
	  return matCreate (side,side,1,cell);
	}
     else return matEmpty(side,side);
   }

