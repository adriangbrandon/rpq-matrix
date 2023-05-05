
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "utilstime.h"

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

// creates an identity matrix

matrix matId (uint64_t side)

{ matrix mat = (matrix)myalloc(sizeof(struct s_matrix));
    uint64_t *data;
    uint64_t nodes;
    side = side ? side : 1;
    mat->width = side;
    mat->height = side;
    mat->logside = numbits(side-1);
    mat->elems = side;
    mat->transposed = 0;
    nodes = (((uint64_t)1) << mat->logside)-1;
    data = (uint64_t*)myalloc(((4*nodes+w-1)/w)*sizeof(uint64_t));
    memset(data,0x99,(nodes+1)/2);
    mat->tree = k2createFrom(mat->logside,4*nodes,data,1);
    return mat;
}

// creates a new copy of A, with its own data

matrix matCopy (matrix A)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    *M = *A;
    if (A->tree != NULL) M->tree = k2copy(A->tree);
    return M;
}

// transposes matrix M

void matTranspose (matrix M)

{ M->transposed = 1-M->transposed;
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

static uint mapId[] = { 0, 1, 2, 3 };
static uint mapTr[] = { 0, 2, 1, 3 };

static uint *mapA,*mapB;

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

typedef struct {
    uint64_t *tree,*levels;
    uint64_t elems,len;
} partition;

static partition empty = { NULL, NULL, 0, 0 };

static partition compose (partition *part, uint level)

{ uint64_t ptr,old;
    uint l,v;
    partition answer;
    uint64_t ptrs[4];

    answer.elems = part[0].elems+part[1].elems+part[2].elems+part[3].elems;
    if (answer.elems == 0) answer = empty;
    else { answer.len = 4 + part[0].len+part[1].len+part[2].len+part[3].len;
        answer.levels = (uint64_t*)myalloc(level*sizeof(uint64_t));
        answer.levels[level-1] = 1;
        answer.tree = (uint64_t*)myalloc((1+(answer.len+w-1)/w)
                                         *sizeof(uint64_t)); //+1 for copyBits
        for (v=0;v<4;v++)
        { bitsWriteA(answer.tree,v,part[v].elems);
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
    for (v=0;v<4;v++)
    { myfree(part[v].tree);
        myfree (part[v].levels);
    }
    return answer;
}

static partition k2sumRC (k2tree treeA,k2node nodeA,k2tree treeB,k2node nodeB,
                          uint level, uint64_t row, uint64_t col, uint has)

{ partition part[4];
    k2node childA[4],childB[4];
    uint valA[4],valB[4];
    uint v,val;
    uint64_t lim;
    partition answer;

    lim = ((uint64_t)1)<<(level-1);
    for (v=0;v<4;v++) valA[v] = valB[v] = 0;
    if ((row == fullSide) || (row < lim))
    { if ((col == fullSide) || (col < lim))
        { if (has & 1) valA[0] = k2hasChild(treeA,nodeA,mapA[0]);
            if (has & 2) valB[0] = k2hasChild(treeB,nodeB,mapB[0]);
        }
        if ((col == fullSide) || (col >= lim))
        { if (has & 1) valA[1] = k2hasChild(treeA,nodeA,mapA[1]);
            if (has & 2) valB[1] = k2hasChild(treeB,nodeB,mapB[1]);
        }
    }
    if ((row == fullSide) || (row >= lim))
    { if ((col == fullSide) || (col < lim))
        { if (has & 1) valA[2] = k2hasChild(treeA,nodeA,mapA[2]);
            if (has & 2) valB[2] = k2hasChild(treeB,nodeB,mapB[2]);
        }
        if ((col == fullSide) || (col >= lim))
        { if (has & 1) valA[3] = k2hasChild(treeA,nodeA,mapA[3]);
            if (has & 2) valB[3] = k2hasChild(treeB,nodeB,mapB[3]);
        }
    }
    // base case, level = 1
    if (level == 1)
    { val =  (valA[0] | valB[0]) +
             ((valA[1] | valB[1]) << 1) +
             ((valA[2] | valB[2]) << 2) +
             ((valA[3] | valB[3]) << 3);
        if (val == 0) return empty;
        answer.elems = popcount(val);
        answer.len = 4;
        answer.tree = (uint64_t*)myalloc(((4+w-1)/w)*sizeof(uint64_t));
        answer.tree[0] = val;
        answer.levels = (uint64_t*)myalloc(1*sizeof(uint64_t));
        answer.levels[0] = 1;
        return answer;
    }

    for (v=0;v<4;v++)
    { if (valA[v]) childA[v] = k2child(treeA,nodeA,mapA[v]);
        if (valB[v]) childB[v] = k2child(treeB,nodeB,mapB[v]);
        part[v] = empty;
    }
    has = valA[0] + 2*valB[0];
    if (has)
        part[0] = k2sumRC (treeA,childA[0],treeB,childB[0],
                           level-1,row,col,has);
    has = valA[1] + 2*valB[1];
    if (has)
        part[1] = k2sumRC (treeA,childA[1],treeB,childB[1],
                           level-1,row,col == fullSide ? col : col-lim, has);
    has = valA[2] + 2*valB[2];
    if (has)
        part[2] = k2sumRC (treeA,childA[2],treeB,childB[2],
                           level-1,row == fullSide ? row : row-lim,col, has);
    has = valA[3] + 2*valB[3];
    if (has)
        part[3] = k2sumRC (treeA,childA[3],treeB,childB[3],
                           level-1,row == fullSide ? row : row-lim,
                           col == fullSide ? col : col-lim, has);

    // combine the 4 quadrants into a single matrix
    return compose(part,level);
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
    M = (matrix)myalloc(sizeof(struct s_matrix));
    M->logside = A->logside;
    matDims(A,NULL,&wA,&hA); matDims(B,NULL,&wB,&hB);
    M->width = mmax(wA,wB); M->height = mmax(hA,hB);
    M->transposed = 0;
    mapA = A->transposed ? mapTr : mapId;
    mapB = B->transposed ? mapTr : mapId;
    sum = k2sum1 (A->tree,B->tree,row,col,&len,&M->elems);
    if (M->elems == 0) M->tree = NULL;
    else M->tree = k2createFrom (k2levels(A->tree),len,sum,1);
    return M;
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
        { part[k].levels = myalloc((level-1)*sizeof(uint64_t));		\
	  part[k].tree = k2merge(part1.tree,part1.len,part2.tree,part2.len,\
			   level-1,&part[k].len,&part[k].elems,part[k].levels);\
	  myfree(part1.tree); myfree(part1.levels);			\
	  myfree(part2.tree); myfree(part2.levels);			\
	}

static partition k2multRC (k2tree treeA, k2node nodeA, k2tree treeB,
                           k2node nodeB, uint level, uint64_t row, uint64_t col)

{
    user_end();
    //fprintf(stdout, "Time %llu\n", (t2-time_t1));
    if(user_diff() > TIMEOUT) {
        //fprintf(stdout, "Time %llu\n", user_diff());
        return empty;
    }

    partition part1,part2,part[4];
    k2node childA[4],childB[4];
    uint valA[4],valB[4];
    uint v,val;
    uint64_t lim;
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
        answer.tree = (uint64_t*)myalloc(((4+w-1)/w)*sizeof(uint64_t));
        answer.tree[0] = val;
        answer.levels = (uint64_t*)myalloc(1*sizeof(uint64_t));
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
    return compose(part,level);
}

static partition k2mult (k2tree treeA, k2node nodeA, k2tree treeB, k2node nodeB,
                         uint level)

{
    user_end();
    //fprintf(stdout, "Time %llu\n", (t2-time_t1));
    if(user_diff() > TIMEOUT) {
        //fprintf(stdout, "Time %llu\n", user_diff());
        return empty;
    }

    partition part1,part2,part[4];
    k2node childA[4],childB[4];
    uint valA[4],valB[4];
    uint v,val;
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
        answer.tree = (uint64_t*)myalloc(((4+w-1)/w)*sizeof(uint64_t));
        answer.tree[0] = val;
        answer.levels = (uint64_t*)myalloc(1*sizeof(uint64_t));
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
    return compose(part,level);
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
        if ((rowA == fullSide) && (colB == fullSide))
            mult = k2mult (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
                           k2levels(A->tree));
        else
            mult = k2multRC (A->tree,k2root(A->tree),B->tree,k2root(B->tree),
                             k2levels(A->tree),rowA,colB);
        myfree (mult.levels);
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
    { Id = matId(mmax(A->width,A->height));
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

{ matrix Id,M,P,S,Ar,E;
    uint64_t elems;
    uint64_t dim = mmax(A->height,A->width);

    if (pos)
    { E = matEmpty(dim,dim);
        Ar = matSum1(A,E,row,fullSide);
        matDestroy(E);
    }
    else
    { Id = matId(dim);
        Ar = matSum1(A,Id,row,fullSide);
        matDestroy(Id);
    }
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

// versions to choose one row or one column, or both

matrix matClos1 (matrix A, uint pos, uint64_t row, uint64_t col)

{ uint64_t nrow,ncol;
    uint64_t side = mmax(A->width,A->height);
    uint64_t cell[2];
    uint64_t test;
    matrix M;
    if (row == fullSide)
    { if (col == fullSide) return matClos(A,pos);
        else
        { matTranspose(A);
            M = matClosRow(A,pos,col,NULL);
            matTranspose(A);
            matTranspose(M);
            return M;
        }
    }
    else
    { if (col == fullSide) return matClosRow(A,pos,row,NULL);
    }
    // both row and col
    nrow = matCollect (A,0,fullSide,col,col,NULL);
    ncol = matCollect (A,row,row,0,fullSide,NULL);
    if (ncol < nrow)
    { test = row;
        matTranspose(A);
        matClosRow(A,pos,col,&test);
        matTranspose(A);
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

