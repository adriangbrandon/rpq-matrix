
// 64-bit version, though matrix sides still must fit in 31 bits

#include <stdlib.h>
#include <string.h>

#include "baseline/matrix.h"
#include "baseline/heap.h"
#include "baseline/hash.h"
#include "baseline/utilstime.h"
#define w32 32

// creates matrix of width x height with n cells (2n ints row,col)
// reorders cells array

// builds row or col

static void sort(uint* r, int lo, int up )

{int i, j, rnd;
    uint tempr;
    while ( up>lo ) {
        i = lo;
        j = up;
        rnd = lo + (random() % (up-lo+1));
        tempr = r[rnd]; r[rnd] = r[lo]; r[lo] = tempr;
        /*** Split file in two ***/
        while ( i<j ) {
            for ( ; r[j] > tempr; j-- );
            for ( r[i]=r[j]; i<j && r[i]<=tempr; i++ );
            r[j] = r[i];
        }
        r[i] = tempr;
        /*** Sort recursively, the smallest first ***/
        if ( i-lo < up-i ) { sort(r,lo,i-1);  lo = i+1; }
        else    { sort(r,i+1,up);  up = i-1; }
    }
}

static uint buildDim (uint size, uint nelems, uint *cells, uint fst, uint snd,
                      uint **ids, uint64_t **pos, uint *other)

{ uint i,d,t;
    uint *dim;

    dim = (uint*)myalloc((1+size)*sizeof(uint));
    for (d=0;d<=size;d++) dim[d] = 0;
    for (i=0;i<nelems;i++) dim[1+cells[2*i+fst]]++;
    t = 0;
    for (d=1;d<=size;d++)
        if (dim[d]) t++;
    *ids = (uint*)myalloc((t+1)*sizeof(uint));
    *pos = (uint64_t*)myalloc((t+1)*sizeof(uint64_t));
    t = 0;
    for (d=1;d<=size;d++)
    { if (dim[d]) (*ids)[t++] = d-1;
        dim[d] += dim[d-1];
    }
    for (i=0;i<nelems;i++) other[dim[cells[2*i+fst]]++] = cells[2*i+snd];
    for (d=0;d<t;d++)
    { if ((*ids)[d] == 0) (*pos)[d] = 0; // could be done a bit better...
        else (*pos)[d] = dim[(*ids)[d]-1];
    }
    (*ids)[t] = size;
    (*pos)[t] = nelems;
    myfree(dim);
    for (d=0;d<t;d++)
    { sort(other+(*pos)[d],0,(*pos)[d+1]-(*pos)[d]-1);
    }
    return t;
}

// packs array into bits bits per element and returns it reallocated

static uint64_t *pack64 (uint64_t *array, uint n, uint bits)

{ uint i;
    for (i=0;i<n;i++)
        write64(array,i,array[i],bits); // should rewrite correctly
    return realloc(array,packedwords(n,bits)*sizeof(uint64_t));
}

// unpacks array and returns unpacked version

static uint64_t *unpack64 (uint64_t *parray, uint n, uint bits)

{ uint i;
    uint64_t *array = (uint64_t*)myalloc(n*sizeof(uint64_t));
    for (i=0;i<n;i++)
        array[i] = access64(parray,i,bits);
    return array;
}

matrix matCreate (uint height, uint width, uint64_t n, uint *cells)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    M->height = height;
    M->width = width;
    M->elems = n;
    if (n == 0)
    { M->rowids = M->colids = NULL;
        M->rowpos = M->colpos = NULL;
        M->rowsbycol = M->colsbyrow = NULL;
        M->nrows = M->ncols = 0;
    }
    else
    { M->colsbyrow = (uint*)myalloc(n*sizeof(uint));
        M->rowsbycol = (uint*)myalloc(n*sizeof(uint));
        M->nrows = buildDim(height,n,cells,0,1,&M->rowids,&M->rowpos,
                            M->colsbyrow);
        M->rowpos = pack64(M->rowpos,M->nrows+1,numbits(M->elems));
        M->ncols = buildDim(width,n,cells,1,0,&M->colids,&M->colpos,
                            M->rowsbycol);
        M->colpos = pack64(M->colpos,M->ncols+1,numbits(M->elems));
    }
    return M;
}


// destroys matrix
void matDestroy (matrix M)

{ if (M == NULL) return;
    myfree (M->rowids); myfree (M->rowpos); myfree (M->colsbyrow);
    myfree (M->colids); myfree (M->colpos); myfree (M->rowsbycol);
    myfree(M);
}

// creates an empty matrix
matrix matEmpty (uint height, uint width)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    M->height = height;
    M->width = width;
    M->elems = 0;
    M->rowids = M->colids = NULL;
    M->rowpos = M->colpos = NULL;
    M->rowsbycol = M->colsbyrow = NULL;
    M->nrows = M->ncols = 0;
    return M;
}

// creates an identity matrix
matrix matId (uint side)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    uint i;
    M->height = M->width = side;
    M->elems = M->nrows = M->ncols = side;
    M->rowids = (uint*)myalloc((side+1)*sizeof(uint));
    M->rowpos = (uint64_t*)myalloc((side+1)*sizeof(uint64_t));
    M->colids = (uint*)myalloc((side+1)*sizeof(uint));
    M->colpos = (uint64_t*)myalloc((side+1)*sizeof(uint64_t));
    M->colsbyrow = (uint*)myalloc(side*sizeof(uint));
    M->rowsbycol = (uint*)myalloc(side*sizeof(uint));
    for (i=0;i<side;i++)
    { M->rowids[i] = M->rowpos[i] = M->colsbyrow[i] = i;
        M->colids[i] = M->colpos[i] = M->rowsbycol[i] = i;
    }
    M->rowids[side] = M->colids[side] = side;
    M->rowpos[side] = M->colpos[side] = side;
    M->rowpos = pack64(M->rowpos,M->nrows+1,numbits(M->elems));
    M->colpos = pack64(M->colpos,M->ncols+1,numbits(M->elems));
    return M;
}

// creates a new copy of M, with its own data
matrix matCopy (matrix M)

{ matrix N = (matrix)myalloc(sizeof(struct s_matrix));
    uint64_t size;
    *N = *M;
    if (N->elems > 0)
    { N->rowids = (uint*)myalloc((N->nrows+1)*sizeof(uint));
        memcpy (N->rowids,M->rowids,(N->nrows+1)*sizeof(uint));
        size = packedwords(M->nrows+1,numbits(M->elems))*sizeof(uint64_t);
        N->rowpos = (uint64_t*)myalloc(size);
        memcpy (N->rowpos,M->rowpos,size);
        N->colsbyrow = (uint*)myalloc(N->elems*sizeof(uint));
        memcpy(N->colsbyrow,M->colsbyrow,N->elems*sizeof(uint));
        N->colids = (uint*)myalloc((N->ncols+1)*sizeof(uint));
        memcpy (N->colids,M->colids,(N->ncols+1)*sizeof(uint));
        size = packedwords(M->ncols+1,numbits(M->elems))*sizeof(uint64_t);
        N->colpos = (uint64_t*)myalloc(size);
        memcpy (N->colpos,M->colpos,size);
        N->rowsbycol = (uint*)myalloc(N->elems*sizeof(uint));
        memcpy(N->rowsbycol,M->rowsbycol,N->elems*sizeof(uint));
    }
    return N;
}

// transpose a matrix, creating a non-allocated copy that shares the
// data. You need not (and should not) matDestroy this copy
// use *M = matTranspose(M) to actually transpose M


struct s_matrix matTranspose (matrix M)

{ uint aux;
    uint *paux;
    struct s_matrix T;
    T.elems = M->elems;
    T.width = M->height;
    T.height = M->width;
    T.nrows = M->ncols;
    T.rowids = M->colids;
    T.rowpos = M->colpos;
    T.colsbyrow = M->rowsbycol;
    T.ncols = M->nrows;
    T.colids = M->rowids;
    T.colpos = M->rowpos;
    T.rowsbycol = M->colsbyrow;
    return T;
}

// writes M to file, which must be opened for writing
void matSave (matrix M, FILE *file)

{ uint64_t num;
    fwrite(&M->elems,sizeof(uint64_t),1,file);
    fwrite(&M->width,sizeof(uint),1,file);
    fwrite(&M->height,sizeof(uint),1,file);
    if (M->elems == 0) return;
    fwrite(&M->nrows,sizeof(uint),1,file);
    fwrite(&M->ncols,sizeof(uint),1,file);

    fwrite(M->rowids,sizeof(uint),M->nrows+1,file);
    num = packedwords(M->nrows+1,numbits(M->elems));
    fwrite(M->rowpos,sizeof(uint64_t),num,file);
    fwrite(M->colsbyrow,sizeof(uint),M->elems,file);

    fwrite(M->colids,sizeof(uint),M->ncols+1,file);
    num = packedwords(M->ncols+1,numbits(M->elems));
    fwrite(M->colpos,sizeof(uint64_t),num,file);
    fwrite(M->rowsbycol,sizeof(uint),M->elems,file);
}

// loads matrix from file, which must be opened for reading
matrix matLoad (FILE *file)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    uint64_t num;
    fread(&M->elems,sizeof(uint64_t),1,file);
    fread(&M->width,sizeof(uint),1,file);
    fread(&M->height,sizeof(uint),1,file);
    if (M->elems == 0)
    { M->rowids = M->colids = NULL;
        M->rowpos = M->colpos = NULL;
        M->rowsbycol = M->colsbyrow = NULL;
        M->nrows = M->ncols = 0;
    }
    else
    { fread(&M->nrows,sizeof(uint),1,file);
        fread(&M->ncols,sizeof(uint),1,file);

        M->rowids = (uint*)myalloc((M->nrows+1)*sizeof(uint));
        fread(M->rowids,sizeof(uint),M->nrows+1,file);
        num = packedwords(M->nrows+1,numbits(M->elems));
        M->rowpos = (uint64_t*)myalloc(num*sizeof(uint64_t));
        fread(M->rowpos,sizeof(uint64_t),num,file);
        M->colsbyrow = (uint*)myalloc(M->elems*sizeof(uint));
        fread(M->colsbyrow,sizeof(uint),M->elems,file);

        M->colids = (uint*)myalloc((M->ncols+1)*sizeof(uint));
        fread(M->colids,sizeof(uint),M->ncols+1,file);
        num = packedwords(M->ncols+1,numbits(M->elems));
        M->colpos = (uint64_t*)myalloc(num*sizeof(uint64_t));
        fread(M->colpos,sizeof(uint64_t),num,file);
        M->rowsbycol = (uint*)myalloc(M->elems*sizeof(uint));
        fread(M->rowsbycol,sizeof(uint),M->elems,file);
    }
    return M;
}

// space of the matrix, in 64-bit words

uint64_t matSpace (matrix M)

{ return (sizeof(struct s_matrix) +
          (M->nrows+1)*sizeof(uint) +
          packedwords(M->nrows+1,numbits(M->elems))*sizeof(uint64_t) +
          M->elems*sizeof(uint) +
          (M->ncols+1)*sizeof(uint) +
          packedwords(M->ncols+1,numbits(M->elems))*sizeof(uint64_t) +
          M->elems*sizeof(uint))
         / (w/8);
}

// dimensions of M, returns #elems and writes the others if not null
uint64_t matDims (matrix M, uint *width, uint *height)

{ if (width) *width = M->width;
    if (height) *height = M->height;
    return M->elems;
}

// accesses a cell

// finds first >= key in r[low..high], may return high+1
static int search(uint key, int low, int high, uint *r)

// gives first >= key
{ int i;
    low--; high++;
    while (high-low > 1)
    {
        i = (high+low) / 2;
        if ( key <= r[i] )  high = i;
        else           low  = i;
    }
    return high;
}

uint matAccess (matrix M, uint row, uint col)

{  // search by row, arbitrarily
    int rpos,rcol;
    uint bits;
    uint64_t from,next;
    uint dif;
    if (M->elems == 0) return 0;
    rpos = search(row,0,M->nrows-1,M->rowids);
    if ((rpos == M->nrows) || (M->rowids[rpos] != row)) return 0;
    bits = numbits(M->elems);
    from = access64(M->rowpos,rpos,bits);
    dif = access64(M->rowpos,rpos+1,bits)-from;
    rcol = search(col,0,dif-1,M->colsbyrow+from);
    if ((rcol == dif) || (M->colsbyrow[from+rcol] != col)) return 0;
    return 1;
}

// recovers all the cells in [r1..r2] x [c1..c2]
// writes 2n integers in buffer, which must be preallocated
// to 2*elems of uint (worst case)
// returns number of elements. just counts if buffer is NULL

uint64_t matCollect (matrix M, uint r1, uint r2,
                     uint c1, uint c2, uint *buffer)

{ int rpos,rcol;
    uint nrows,rowid,erow;
    uint64_t srow;
    uint dif,bits;
    uint64_t t = 0;
    if (M->elems == 0) return 0;
    bits = numbits(M->elems);
    nrows = M->nrows;
    rpos = search(r1,0,nrows-1,M->rowids);
    while ((rpos < nrows) && ((rowid = M->rowids[rpos]) <= r2))
    { srow = access64(M->rowpos,rpos,bits);
        dif = access64(M->rowpos,rpos+1,bits)-srow;
        rcol = search(c1,0,dif-1,M->colsbyrow+srow);
        if (buffer)
        { while ((rcol < dif) && (M->colsbyrow[srow+rcol] <= c2))
            { buffer[2*t] = rowid;
                buffer[2*t+1] = M->colsbyrow[srow+rcol];
                t++;
                rcol++;
            }
        }
        else if (rcol < dif)
        { t -= rcol;
            rcol = search(c2+1,0,dif-1,M->colsbyrow+srow);
            t += rcol;
        }
        rpos++;
    }
    return t;
}

// (boolean) sum of two matrices, assumed to be of the same side

static uint merge (uint *dataA, uint fa, uint *dataB, uint fb, uint *dataM)

{ uint c,ta,tb;
    c = ta = tb = 0;
    while ((ta < fa) && (tb < fb))
    { if (dataA[ta] < dataB[tb]) dataM[c++] = dataA[ta++];
        else if (dataA[ta] > dataB[tb]) dataM[c++] = dataB[tb++];
        else
        { dataM[c++] = dataA[ta];
            ta++; tb++;
        }
    }
    if (ta < fa)
    { memcpy(dataM+c,dataA+ta,(fa-ta)*sizeof(uint));
        c += fa-ta;
    }
    else if (tb < fb)
    { memcpy(dataM+c,dataB+tb,(fb-tb)*sizeof(uint));
        c += fb-tb;
    }
    return c;
}

static uint matMerge (uint na, uint nb, uint64_t *posA, uint64_t *posB,
                      uint64_t *posM, uint *idsA, uint *idsB, uint *idsM,
                      uint *dataA, uint *dataB, uint *dataM)

{ uint p,pa,pb,dif;
    uint64_t c;
    c = 0;
    p = pa = pb = 0;
    while ((pa < na) || (pb < nb))
    { posM[p] = c;
        if (idsA[pa] < idsB[pb])
        { idsM[p] = idsA[pa];
            dif = posA[pa+1]-posA[pa];
            memcpy(dataM+c,dataA+posA[pa],dif*sizeof(uint));
            pa++; c += dif;
        }
        else if (idsB[pb] < idsA[pa])
        { idsM[p] = idsB[pb];
            dif = posB[pb+1]-posB[pb];
            memcpy(dataM+c,dataB+posB[pb],dif*sizeof(uint));
            pb++; c += dif;
        }
        else
        { idsM[p] = idsA[pa];
            c += merge(dataA+posA[pa],posA[pa+1]-posA[pa],
                       dataB+posB[pb],posB[pb+1]-posB[pb],dataM+c);
            pa++; pb++;
        }
        p++;
    }
    posM[p] = c;
    return p;
}

matrix matSum (matrix A, matrix B)

{ matrix M;
    uint p,pa,pb,na,nb,c;
    uint64_t *posA,*posB;
    uint bitsA,bitsB,bitsM;

    if ((A->elems == 0) || (B->elems == 0))
    { if (A->elems == 0) M = matCopy(B); else M = matCopy(A);
        M->width = mmax(A->width,B->width);
        M->height = mmax(A->height,B->height);
        return M;
    }
    M = (matrix)myalloc(sizeof(struct s_matrix));
    M->width = mmax(A->width,B->width);
    M->height = mmax(A->height,B->height);
    bitsA = numbits(A->elems);
    bitsB = numbits(B->elems);

    M->rowids = (uint*)myalloc((A->nrows+B->nrows+1)*sizeof(uint));
    M->rowpos = (uint64_t*)myalloc((A->nrows+B->nrows+1)*sizeof(uint64_t));
    M->colsbyrow = (uint*)myalloc((A->elems+B->elems)*sizeof(uint));
    posA = unpack64(A->rowpos,A->nrows+1,bitsA);
    posB = unpack64(B->rowpos,B->nrows+1,bitsB);
    M->nrows = matMerge (A->nrows,B->nrows,posA,posB,M->rowpos,
                         A->rowids,B->rowids,M->rowids,
                         A->colsbyrow,B->colsbyrow,M->colsbyrow);
    myfree(posA); myfree(posB);
    M->rowids[M->nrows] = M->height;
    M->rowids = (uint*)myrealloc(M->rowids,(M->nrows+1)*sizeof(uint));
    M->elems = M->rowpos[M->nrows];
    bitsM = numbits(M->elems);
    M->rowpos = pack64(M->rowpos,M->nrows+1,bitsM);
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,M->elems*sizeof(uint));

    M->colids = (uint*)myalloc((A->ncols+B->ncols+1)*sizeof(uint));
    M->colpos = (uint64_t*)myalloc((A->ncols+B->ncols+1)*sizeof(uint64_t));
    posA = unpack64(A->colpos,A->ncols+1,bitsA);
    posB = unpack64(B->colpos,B->ncols+1,bitsB);
    M->rowsbycol = (uint*)myalloc(M->elems*sizeof(uint));
    M->ncols = matMerge (A->ncols,B->ncols,posA,posB,M->colpos,
                         A->colids,B->colids,M->colids,
                         A->rowsbycol,B->rowsbycol,M->rowsbycol);
    myfree(posA); myfree(posB);
    M->colids[M->ncols] = M->width;
    M->colids = (uint*)myrealloc(M->colids,(M->ncols+1)*sizeof(uint));
    M->colpos = pack64(M->colpos,M->ncols+1,bitsM);

    return M;
}

// version with one row or one column, or both

static matrix matOne (uint height, uint width, uint row, uint col)

{ matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    uint bits;
    M->elems = M->nrows = M->ncols = 1;
    M->height = height; M->width = width;
    bits = numbits(M->elems);

    M->rowids = (uint*)myalloc(2*sizeof(uint));
    M->rowids[0] = row; M->rowids[1] = height;
    M->rowpos = (uint64_t*)myalloc(packedwords(2,bits)*sizeof(uint64_t));
    write64(M->rowpos,0,0,bits); write64(M->rowpos,1,1,bits);
    M->colsbyrow = (uint*)myalloc(1*sizeof(uint));
    M->colsbyrow[0] = col;

    M->colids = (uint*)myalloc(2*sizeof(uint));
    M->colids[0] = col; M->colids[1] = width;
    M->colpos = (uint64_t*)myalloc(packedwords(2,bits)*sizeof(uint64_t));
    write64(M->colpos,0,0,bits); write64(M->colpos,1,1,bits);
    M->rowsbycol = (uint*)myalloc(1*sizeof(uint));
    M->rowsbycol[0] = row;

    return M;
}

static void copy1 (uint **colids, uint64_t **colpos, uint **rowsbycol,
                   uint **rowids, uint64_t **rowpos, uint **colsbyrow,
                   uint *data, uint col, uint width, uint height, uint nelems)

{ uint i;
    uint bits = numbits(nelems);
    *colids = (uint*)myalloc(2*sizeof(uint));
    *colpos = (uint64_t*)myalloc(packedwords(2,bits)*sizeof(uint64_t));
    (*colids)[0] = col; (*colids)[1] = width;
    write64(*colpos,0,0,bits); write64(*colpos,1,nelems,bits);
    *rowsbycol = (uint*)myalloc(nelems*sizeof(uint));
    memcpy(*rowsbycol,data,nelems*sizeof(uint));
    *rowids = (uint*)myalloc((nelems+1)*sizeof(uint));
    *rowpos = (uint64_t*)myalloc(packedwords(nelems+1,bits)*sizeof(uint64_t));
    memcpy(*rowids,data,nelems*sizeof(uint));
    (*rowids)[nelems] = height;
    *colsbyrow = (uint*)myalloc(nelems*sizeof(uint));
    for (i=0;i<=nelems;i++) write64(*rowpos,i,i,bits);
    for (i=0;i<nelems;i++) (*colsbyrow)[i] = col;
}

matrix matSum1 (uint row, matrix A, matrix B, uint col)

{ uint posA,posB,i;
    uint height = mmax(A->height,B->height);
    uint width = mmax(A->width,B->width);
    matrix M;
    struct s_matrix At,Bt;
    uint bitsA,bitsB;
    uint64_t pa,pb;
    uint na,nb;
    uint *data;

    if (row == fullSide)
        if (col == fullSide)
            return matSum(A,B);
        else
        { posA = search(col,0,A->ncols-1,A->colids);
            posB = search(col,0,B->ncols-1,B->colids);
            if (((posA >= A->ncols) || (A->colids[posA] != col)) &&
                ((posB >= B->ncols) || (B->colids[posB] != col)))
                return matEmpty (height,width);
            M = (matrix)myalloc(sizeof(struct s_matrix));
            M->height = height; M->width = width;
            M->ncols = 1;
            if ((posA >= A->ncols) || (A->colids[posA] != col))
            { bitsB = numbits(B->elems);
                M->elems = access64(B->colpos,posB+1,bitsB) -
                           access64(B->colpos,posB,bitsB);
                copy1 (&M->colids,&M->colpos,&M->rowsbycol,
                       &M->rowids,&M->rowpos,&M->colsbyrow,
                       B->rowsbycol+access64(B->colpos,posB,bitsB),
                       col,width,height,M->elems);
            }
            else if ((posB >= B->ncols) || (B->colids[posB] != col))
            { bitsA = numbits(A->elems);
                M->elems = access64(A->colpos,posA+1,bitsA) -
                           access64(A->colpos,posA,bitsA);
                copy1 (&M->colids,&M->colpos,&M->rowsbycol,
                       &M->rowids,&M->rowpos,&M->colsbyrow,
                       A->rowsbycol+access64(A->colpos,posA,bitsA),
                       col,width,height,M->elems);
            }
            else
            { bitsA = numbits(A->elems);
                pa = access64(A->colpos,posA,bitsA);
                na = access64(A->colpos,posA+1,bitsA)-pa;
                bitsB = numbits(B->elems);
                pb = access64(B->colpos,posB,bitsB);
                nb = access64(B->colpos,posB+1,bitsB)-pb;
                data = (uint*)myalloc((na+nb)*sizeof(uint));
                M->elems = merge(A->rowsbycol+pa,na,
                                 B->rowsbycol+pb,nb,data);
                copy1 (&M->colids,&M->colpos,&M->rowsbycol,
                       &M->rowids,&M->rowpos,&M->colsbyrow,
                       data,col,width,height,M->elems);
                myfree(data);
            }
            M->nrows = M->elems;
            return M;
        }
    else if (col == fullSide) // sorry for this, I couldn't resist
    { At = matTranspose(A); Bt = matTranspose(B);
        M = matSum1(col,&At,&Bt,row);
        *M = matTranspose(M);
        return M;
    }
    else // both row and col restrictions
    { if (matAccess(A,row,col) && matAccess(B,row,col))
            return matOne (height,width,row,col);
        else return matEmpty (height,width);
    }
}

// (boolean) product of two matrices, assumed to be of the same side
// only row of A and col of B are considered if not fullSide

// assumes na << nb and data is preallocated (presumably to na elems)

static uint inter120 (uint *dataA, uint na, uint *dataB, uint nb)

{ int i,j;
    // uint p = 0;
    j = 0;
    for (i=0;i<na;i++)
    { uint d = dataA[i];
        uint s = 1;
        j--;
        while ((j+s < nb) && (dataB[j+s] < d)) s <<= 1;
        j = search(d,j+1+(s>>1),mmin(j+s,nb-1),dataB);
        if ((j < nb) && (dataB[j] == d)) return 1; // data[p++] = d;
    }
    return 0; // p;
}

static uint inter0 (uint *dataA, uint na, uint *dataB, uint nb)

{ int i,j;
    // uint p = 0;
    i = j = 0;
    while ((i < na) && (j < nb))
    { uint da = dataA[i];
        uint db = dataB[j];
        if (da == db) return 1; // data[p++] = da;
        if (da < db) i++; else j++;
    }
    return 0; // p;
}

static uint inters0 (uint *dataA, uint na, uint *dataB, uint nb)

{ if (nb > 30*na) return inter120(dataA,na,dataB,nb);
    if (na > 30*nb) return inter120(dataB,nb,dataA,na);
    return inter0(dataA,na,dataB,nb);
}

typedef struct
{ uint rowB,colA;
    uint ncols,nrows;
    uint *cbyr,*rbyc;
    uint pcbyr;
} ttask;


// assumes task is preallocated, to at least min(na,nb) elements

static uint inter12 (uint *dataA, uint na, uint *dataB, uint nb,
                     ttask *task, uint aisr)

{ int i,j;
    uint p = 0;
    j = 0;
    for (i=0;i<na;i++)
    { uint d = dataA[i];
        uint s = 1;
        j--;
        while ((j+s < nb) && (dataB[j+s] < d)) s <<= 1;
        j = search(d,j+1+(s>>1),mmin(j+s,nb-1),dataB);
        if ((j < nb) && (dataB[j] == d))
        { if (aisr) { task[p].colA = i; task[p].rowB = j; }
            else { task[p].colA = j; task[p].rowB = i; }
            p++; j++;
        }
    }
    return p;
}

static uint inter (uint *dataA, uint na, uint *dataB, uint nb, ttask *task)

{ int i,j;
    uint p = 0;
    i = j = 0;
    while ((i < na) && (j < nb))
    { uint da = dataA[i];
        uint db = dataB[j];
        if (da == db)
        { task[p].colA = i; task[p].rowB = j;
            p++; i++; j++;
        }
        else if (da < db) i++;
        else j++;
    }
    return p;
}

static uint inters (uint *dataA, uint na, uint *dataB, uint nb, ttask *task)

{ if (nb > 30*na) return inter12(dataA,na,dataB,nb,task,1);
    if (na > 30*nb) return inter12(dataB,nb,dataA,na,task,0);
    return inter(dataA,na,dataB,nb,task);
}



matrix matMult (matrix A, matrix B)

{
#if TIMEOUT
    if(time_diff() > TIMEOUT) {
        //fprintf(stdout, "Time %llu\n", user_diff());
        return matEmpty(A->height, B->width);
    }
#endif

    ttask *task;
    uint c,i,j,m,p,pr,pc;
    heap Hr,Hc;
    uint size;
    uint *colc; // column counter to later build the transpose
    uint64_t pos,prc;
    uint bitsA,bitsB;
    uint *stackp,*stackv; // implement colc as an initializable array
    uint stacks;

    // tasks are columns of A that are also rows of B
    task = (ttask*)myalloc(mmin(A->ncols,B->nrows)*sizeof(ttask));
    p = inters(A->colids,A->ncols,B->rowids,B->nrows,task);
    if (p == 0) // empty result!
    { myfree(task);
        return matEmpty(A->height,B->width);
    }
    task = (ttask*)myrealloc(task,p*sizeof(ttask));
    Hr = (heap)myalloc(p*sizeof(struct s_heap)); // to merge rows
    Hc = (heap)myalloc(p*sizeof(struct s_heap)); // to merge columns

    // we first fill M in row-wise form
    matrix M = (matrix)myalloc(sizeof(struct s_matrix));
    M->rowids = (uint*)myalloc((A->nrows+1)*sizeof(uint));
    M->rowpos = (uint64_t*)myalloc((A->nrows+1)*sizeof(uint64_t));
    size = 1 << numbits(A->nrows);
    M->colsbyrow = (uint*)myalloc(size*sizeof(uint)); // growing array
    pr = prc = 0;
    colc = (uint*)myalloc(B->width*sizeof(uint));
    stackp = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stackv = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stacks = 0;
    // for (c=0;c<B->width;c++) colc[c] = 0;

    // create the p tasks
    bitsA = numbits(A->elems);
    bitsB = numbits(B->elems);
    for (i=0;i<p;i++)
    { pos = access64(B->rowpos,task[i].rowB,bitsB);
        task[i].ncols = access64(B->rowpos,task[i].rowB+1,bitsB)-pos;
        task[i].cbyr = B->colsbyrow+pos;
        pos = access64(A->colpos,task[i].colA,bitsA);
        task[i].nrows = access64(A->colpos,task[i].colA+1,bitsA)-pos;
        task[i].rbyc = A->rowsbycol+pos;
        Hr[i].key = *task[i].rbyc++;
        task[i].nrows--;
        Hr[i].data = i;
    }

    // traverse them by increasing row
    heapify(Hr,p);
    while (p)
    { struct s_heap hr;
        // get all equal mminimum hr.key and collect hr.data in Hc
        m = 0;
        do { hr = findMin(Hr,p);
            i = hr.data;
            // store for column merge later
            Hc[m].key = *task[i].cbyr; // store for column merge
            Hc[m].data = i;
            task[i].pcbyr = 1; // first is already used
            m++;
            // replace by the next row
            if (task[i].nrows)
            { replaceMin(Hr,p,*task[i].rbyc++,i);
                task[i].nrows--;
            }
            else extractMin(Hr,p--);
        }
        while (p && (findMin(Hr,p).key == hr.key));
        // create new row in M and merge its columns
        M->rowids[pr] = hr.key;
        M->rowpos[pr] = prc;
        pr++;
        heapify(Hc,m);
        while (m)
        { struct s_heap hc = findMin(Hc,m);
            uint k = hc.key;
            i = hc.data;
            if (prc == size) // realloc colsbyrow
            { size <<= 1;
                M->colsbyrow = (uint*)myrealloc(M->colsbyrow,
                                                size*sizeof(uint));
            }
            M->colsbyrow[prc++] = k;
            if ((colc[k] >= stacks) || (stackp[colc[k]] != k)) // invalid value
            {   stackp[stacks] = k;
                stackv[stacks] = 1;
                colc[k] = stacks++;
            }
            else stackv[colc[k]]++;
            do {   // replace by the next col
                if (task[i].pcbyr < task[i].ncols)
                    replaceMin(Hc,m,task[i].cbyr[task[i].pcbyr++],i);
                else { extractMin(Hc,m--);
                    if (!m) break;
                }
                hc = findMin(Hc,m);
                i = hc.data;
            }
            while (hc.key == k);
        }
    }

    // close this stage
    myfree(task);
    myfree(Hr); myfree(Hc);
    if (prc == 0) // empty result!
    { myfree(M->rowids); myfree(M->rowpos); myfree(M->colsbyrow); myfree(M);
        myfree(colc);
        myfree(stackv); myfree(stackp);
        return matEmpty(A->height,B->width);
    }
    M->elems = prc;
    M->nrows = pr;
    M->rowids[pr] = A->height;
    M->rowpos[pr] = prc;
    M->rowids = (uint*)myrealloc(M->rowids,(pr+1)*sizeof(uint));
    M->rowpos = (uint64_t*)myrealloc(M->rowpos,(pr+1)*sizeof(uint64_t));
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,prc*sizeof(uint));

    // now create the transposed representation

    //M->colids = (uint*)myalloc((B->width+1)*sizeof(uint));
    //pc = 0;
    //for (c=0;c<B->width;c++)
    //{ if (colc[c]) M->colids[pc++] = c; }
    M->colids = stackp;
    for (pc=0;pc<stacks;pc++){
        colc[stackp[pc]] = stackv[pc];
    }
    sort(M->colids,0,pc-1);
    M->colids[pc] = B->width;
    M->ncols = pc;
    M->colids = (uint*)myrealloc(M->colids,(pc+1)*sizeof(uint));
    M->colpos = (uint64_t*)myalloc((pc+1)*sizeof(uint64_t));
    M->rowsbycol = (uint*)myalloc(prc*sizeof(uint));
    M->colpos[0] = 0;
    for (c=1;c<pc;c++)
        M->colpos[c] = M->colpos[c-1] + colc[M->colids[c-1]];
    M->colpos[pc] = prc;
    for (c=0;c<pc;c++) colc[M->colids[c]] = M->colpos[c];

    for (pr=0;pr<M->nrows;pr++)
    { uint f = M->rowpos[pr+1];
        uint r = M->rowids[pr];
        for (prc=M->rowpos[pr];prc<f;prc++)
            M->rowsbycol[colc[M->colsbyrow[prc]]++] = r;
    }
    myfree(colc);
    myfree(stackv);
    M->rowpos = pack64(M->rowpos,pr+1,numbits(M->elems));
    M->colpos = pack64(M->colpos,pc+1,numbits(M->elems));
    M->width = B->width;
    M->height = A->height;
    return M;
}

matrix matMult1R (uint r, matrix A, matrix B)

{


    ttask *task;
    uint i,j,p,pr,pc;
    heap Hc;
    matrix M;
    uint bitsA,bitsB,bits;
    uint64_t posA,posB,prc;

    pr = search(r,0,A->nrows,A->rowids);
    if (A->rowids[pr] != r) return matEmpty(A->height,B->width);
    bitsA = numbits(A->elems);
    bitsB = numbits(B->elems);

    // tasks are columns of A that are also rows of B
    task = (ttask*)myalloc(mmin(A->ncols,B->nrows)*sizeof(ttask));
    posA = access64(A->rowpos,pr,bitsA);
    p = inters(A->colsbyrow+posA,access64(A->rowpos,pr+1,bitsA)-posA,
               B->rowids,B->nrows,task);
    if (p == 0) // empty result!
    { myfree(task);
        return matEmpty(A->height,B->width);
    }
    task = (ttask*)myrealloc(task,p*sizeof(ttask));
    Hc = (heap)myalloc(p*sizeof(struct s_heap)); // to merge columns

    M = (matrix)myalloc(sizeof(struct s_matrix));
    M->colsbyrow = (uint*)myalloc(B->ncols*sizeof(uint));
    prc = 0;

    // create the p tasks
    for (i=0;i<p;i++)
    { posB = access64(B->rowpos,task[i].rowB,bitsB);
        task[i].ncols = access64(B->rowpos,task[i].rowB+1,bitsB)-posB;
        task[i].cbyr = B->colsbyrow+posB;
        Hc[i].key = *task[i].cbyr++;
        task[i].ncols--;
        Hc[i].data = i;
    }

    // traverse them by increasing column
    heapify(Hc,p);
    while (p)
    { struct s_heap hc = findMin(Hc,p);
        uint k = hc.key;
        i = hc.data;
        M->colsbyrow[prc++] = k;
        do {   // replace by the next col
            if (task[i].ncols)
            { replaceMin(Hc,p,*task[i].cbyr++,i);
                task[i].ncols--;
            }
            else
            { extractMin(Hc,p--);
                if (!p) break;
            }
            hc = findMin(Hc,p);
            i = hc.data;
        }
        while (hc.key == k);
    }

    // close this stage
    myfree(task);
    myfree(Hc);
    if (prc == 0) // empty result!
    { myfree(M->colsbyrow); myfree(M);
        return matEmpty(A->height,B->width);
    }
    bits = numbits(prc);
    M->nrows = 1;
    M->rowids = (uint*)myalloc(2*sizeof(uint));
    M->rowids[0] = r; M->rowids[1] = A->height;
    M->rowpos = (uint64_t*)myalloc(packedwords(2,bits)*sizeof(uint64_t));
    write64(M->rowpos,0,0,bits); write64(M->rowpos,1,prc,bits);
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,prc*sizeof(uint));

    // create the transposed version
    M->ncols = prc;
    M->colids = (uint*)myalloc((prc+1)*sizeof(uint));
    M->colpos = (uint64_t*)myalloc(packedwords(prc+1,bits)*sizeof(uint64_t));
    M->rowsbycol = (uint*)myalloc(prc*sizeof(uint));
    memcpy(M->colids,M->colsbyrow,prc*sizeof(uint));
    for (pc=0;pc<=prc;pc++) write64(M->colpos,pc,pc,bits);
    for (pc=0;pc<prc;pc++) M->rowsbycol[pc] = r;
    M->colids[pc] = B->width;

    M->width = B->width;
    M->height = A->height;
    M->elems = prc;
    return M;
}

// version with one row or one column, or both

matrix matMult1 (uint row, matrix A, matrix B, uint col)

{ matrix M;
    struct s_matrix At,Bt;
    uint bitsA,bitsB;
    uint64_t posA,posB;

    if ((A->elems == 0) || (B->elems == 0))
        return matEmpty(A->height,B->width);
    if ((row == fullSide) && (col == fullSide)) return matMult(A,B);
    if ((row != fullSide) && (col != fullSide))
    { bitsA = numbits(A->elems); bitsB = numbits(B->elems);
        posA = access64(A->rowpos,row,bitsA);
        posB = access64(B->colpos,col,bitsB);
        if (inters0(A->colsbyrow+posA,access64(A->rowpos,row+1,bitsA)-posA,
                    B->rowsbycol+posB,access64(B->colpos,col+1,bitsB)-posB))
        return matOne(A->height,B->width,row,col);
        else return matEmpty(A->height,B->width);
    }
    if (row == fullSide)
    { At = matTranspose(A); Bt = matTranspose(B);
        M = matMult1R(col,&Bt,&At);
        *M = matTranspose(M);
        return M;
    }
    return matMult1R(row,A,B);
}

// transitive closure of a matrix, pos says if it's + rather than *

matrix matClos0 (matrix A, uint pos)

{ matrix M,P,S,Id;
    uint elems;

    elems = A->elems;
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
    if (!pos)
    { Id = matId(mmax(A->width,A->height));
        M = S; S = matSum(S,Id); matDestroy(M);
        matDestroy(Id);
    }
    return S;
}

// versions to chose row row or column col at the end
// here we restrict first to the row/column even if we lose the
// exponential increase in the paths

// if coltest not null, checks that the cell *coltest exists and
// returns there 1 or 0 as soon as it can
// if not null, I implies starting with I, I x A* or I x A+

static matrix matClosRow (uint row, matrix I, matrix A, uint pos, uint *coltest)

{ matrix M,P,S,E;
    uint64_t elems;
    uint dim = mmax(A->height,A->width);

    if (I == NULL)
    { if (pos) E = matEmpty(dim,dim);
        else E = matOne(dim,dim,row,row);
        S = matSum1 (row,A,E,fullSide);
        matDestroy(E);
    }
    else
    { if (pos) S = matMult1(row,I,A,fullSide);
        else { E = matEmpty(dim,dim);
            S = matSum1(row,I,E,fullSide);
            matDestroy(E);
        }
    }
    elems = S->elems;
    if (coltest && matAccess(S,row,*coltest))
    { matDestroy(S); *coltest = 1; return NULL; }
    P = matMult (S,A);
    M = S; S = matSum (P,S); matDestroy(M);
    while (S->elems != elems)
    { if (coltest && matAccess(S,row,*coltest))
        { matDestroy(S); matDestroy(P); *coltest = 1; return NULL; }
        elems = S->elems;
        M = P; P = matMult(P,A); matDestroy(M);
        M = S; S = matSum(S,P); matDestroy(M);
    }
    matDestroy(P);
    if (coltest) { matDestroy(S); *coltest = 0; return NULL; }
    return S;
}

// versions to choose one row or one column, or both

matrix matClos1 (uint row, matrix A, uint pos, uint col)

{ uint nrow,ncol;
    uint side = mmax(A->width,A->height);
    uint test;
    matrix M;
    struct s_matrix At;

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
    if (test) return matOne (side,side,row,col);
    else return matEmpty(side,side);
}

// computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)

matrix matMultClos1 (uint row, matrix A, matrix B, uint pos, uint col)

{ uint nrow,ncol;
    uint side = mmax(A->width,A->height);
    uint test;
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

matrix matClosMult1 (uint row, matrix A, uint pos, matrix B, uint col)

{ matrix M;
    struct s_matrix At,Bt;
    At = matTranspose(A);
    Bt = matTranspose(B);
    M = matMultClos1(col,&Bt,&At,pos,row);
    *M = matTranspose(M);
    return M;
}


// finds the strongly connected components of the matrix

typedef struct
{ uint index;
    uint lowlink;
    uint compid; // stores onstack in highest bit
} tarjan;

static uint *stack; // stores pointers to rowids
static uint pstack;
static uint ncomp; // components defined up to now
static tarjan *info; // for the algorithm
static uint tindex; // global index
static hashTable hash; // to map rowids to rows

static uint *comptr; // points to the start of each component in compnodes
static uint *compnodes; // positions in info of the nodes grouped by component
static uint pcomp; // pointer to next free pos in compnodes

#define onstack (((uint)1)<<(w32-1)) // bitset for onstack

uint bitsM;

static void strongconnect (matrix M, uint v)

{ uint r,u;
    uint64_t c,next;
    info[v].index  = tindex;
    info[v].lowlink = tindex;
    info[v].compid = onstack;
    tindex++;
    stack[pstack++] = v;
    next = access64(M->rowpos,v+1,bitsM);

    for (c=access64(M->rowpos,v,bitsM); c<next; c++)
    { r = M->colsbyrow[c];
        u = hashSearch(hash,r);
        if (u != noval) // sink nodes are ignored for now
        { if (info[u].index == 0)
            { strongconnect (M,u);
                info[v].lowlink = mmin(info[v].lowlink,info[u].lowlink);
            }
            else if (info[u].compid & onstack)
            { info[v].lowlink = mmin(info[v].lowlink,info[u].index);
            }
        }
    }

    if (info[v].lowlink == info[v].index)
    { comptr[ncomp] = pcomp;
        while (1)
        { u = stack[--pstack];
            info[u].compid = ncomp; // onstack = 0
            compnodes[pcomp++] = M->rowids[u];
            if (u == v) break;
        }
        ncomp++;
    }
}

// max-heap of pairs

typedef struct {
    uint r,c;
} tpair;

#define less(x,y) ((x.r < y.r) || ((x.r == y.r) && (x.c < y.c)))

#define swap(x,y) \
   { uint aux = pairs[2*(x)]; pairs[2*(x)] = pairs[2*(y)]; \
     pairs[2*(y)] = swap; swap = pairs[2*(x)+1]; \
     pairs[2*(x)+1] = pairs[2*(y)+1]; pairs[2*(y)+1] = swap; }

// min-heap of pairs that runs backwards
// n is the heap size, root is H[n-1]
// but indices are given as if it were forwards H[1..n]

static inline void Psiftdown (tpair *H, uint n, uint i)

{ tpair val = H[n-i];
    uint pos;
    while (2*i <= n)
    { if ((2*i+1 > n) || less(H[n-2*i],H[n-(2*i+1)])) pos = 2*i;
        else pos = 2*i+1;
        if (less(H[n-pos],val)) { H[n-i] = H[n-pos]; i = pos; }
        else break;
    }
    H[n-i] = val;
}

// heapsorts edges to get them unique and in order and at the same
// time fills the components of the resulting graph of components
// returns number of rows in such matrix

static uint createReduced (tpair *H, uint64_t n,
                           uint *colsbyrow, uint64_t *rowpos)

{ uint64_t i;
    uint64_t p,b; // writes at colsbyrow[p], heap is H[b..n-1]
    uint rp; // writes at rowpos[rp]
    tpair val;
    // heapify
    for (i=n/2;i>=1;i--) Psiftdown(H,n,i);
    // extracts
    p = b = 0; val = H[n-1];
    rowpos[0] = 0;
    colsbyrow[p++] = val.c;
    rp = 1;
    while (b < n-1)
    { 	// extract min
        H[n-1] = H[b++];
        Psiftdown(H+b,n-b,1);
        if ((H[n-1].r != val.r) || (H[n-1].c != val.c))
        { if (val.r != H[n-1].r)
            { val.r = H[n-1].r;
                rowpos[rp] = p;
                rp++;
            }
            val.c = H[n-1].c;
            colsbyrow[p++] = val.c;
        }
    }
    rowpos[rp] = p;
}

// min-heap of uints that runs backwards

static inline void siftdown (uint *H, uint64_t n, uint64_t i)

{ uint val = H[n-i];
    uint64_t pos;
    while (2*i <= n)
    { if ((2*i+1 > n) || (H[n-2*i] < H[n-(2*i+1)])) pos = 2*i;
        else pos = 2*i+1;
        if (H[n-pos] < val) { H[n-i] = H[n-pos]; i = pos; }
        else break;
    }
    H[n-i] = val;
}

static uint64_t myheapsort (uint *H, uint64_t n)

{ uint64_t i;
    uint64_t p,b; // writes at p, heap is H[b..n-1]
    uint val;
    // heapify
    for (i=n/2;i>=1;i--) siftdown(H,n,i);
    // extracts
    p = b = 0; val = H[n-1];
    while (b < n-1)
    {	// extract min
        H[n-1] = H[b++];
        siftdown(H+b,n-b,1);
        if (H[n-1] != val) { H[p++] = val; val = H[n-1]; }
    }
    H[p++] = val;
    return p;
}

// C is already acyclic, propagate reachability

static uint *colsbyrow;
static uint64_t *rowpos;

#define resize(array,size,bits,newsize) \
   { if ((newsize) > size) \
      { if ((newsize) >> bits) \
	   { bits = numbits(newsize); \
             array=(uint*)myrealloc(array,(((uint64_t)1)<<bits)*sizeof(uint)); \
	   } \
	size = (newsize); \
      } \
   }

static void propagate (struct s_matrix C)

{ uint64_t j,size,pos,lsize;
    uint bits;
    uint v,u;

    size = C.elems;
    bits = numbits(size);
    colsbyrow = (uint*)myalloc((((uint64_t)1)<<bits)*sizeof(uint));
    rowpos = (uint64_t*)myalloc((C.nrows+1)*sizeof(uint64_t));
    pos = 0;

    for (v=0;v<C.nrows;v++)
    { rowpos[v] = pos;
        lsize = C.rowpos[v+1]-C.rowpos[v];
        resize(colsbyrow,size,bits,pos+lsize);
        memcpy(colsbyrow+pos,C.colsbyrow+C.rowpos[v],lsize*sizeof(uint));
        pos += lsize;
        for (j=C.rowpos[v];j<C.rowpos[v+1];j++)
        { u = C.colsbyrow[j];
            if (u < v) // not a sink (u < ncomp) and not v itself
            { lsize = rowpos[u+1]-rowpos[u];
                resize(colsbyrow,size,bits,pos+lsize);
                memcpy(colsbyrow+pos,colsbyrow+rowpos[u],
                       lsize*sizeof(uint));
                pos += lsize;
            }
        }
        // sort list just created and remove repeated elements
        pos = rowpos[v] + myheapsort(colsbyrow+rowpos[v],pos-rowpos[v]);
    }
    rowpos[C.nrows] = pos;
    colsbyrow = (uint*)myrealloc(colsbyrow,pos*sizeof(uint));
}

// expand representatives of components to all their nodes

static uint *expand (uint nrows, uint64_t *nedges)

{ uint *edges;
    uint usize,vsize;
    uint64_t k,size,pos;
    uint v,u,i,j;
    uint bits;

    size = 2*rowpos[nrows];
    bits = numbits(size);
    edges = myalloc((((uint64_t)1)<<bits)*sizeof(uint));
    pos = 0;
    for (v=0;v<nrows;v++)
    {  // add product of the connected components
        vsize = comptr[v+1]-comptr[v];
        for (k=rowpos[v];k<rowpos[v+1];k++)
        { u = colsbyrow[k];
            if (u < ncomp)
            { usize = comptr[u+1]-comptr[u];
                resize(edges,size,bits,pos+2*vsize*(uint64_t)usize);
                for (i=comptr[v];i<comptr[v+1];i++)
                    for (j=comptr[u];j<comptr[u+1];j++)
                    { edges[pos++] = compnodes[i];
                        edges[pos++] = compnodes[j];
                    }
            }
            else // a sink, real id is u-ncomp
            { resize(edges,size,bits,pos+2*(uint64_t)vsize);
                for (i=comptr[v];i<comptr[v+1];i++)
                { edges[pos++] = compnodes[i];
                    edges[pos++] = u-ncomp;
                }
            }
        }
    }
    myfree(rowpos);
    myfree(colsbyrow);
    *nedges = pos/2;
    edges = (uint*)myrealloc(edges,pos*sizeof(uint));
    return edges;
}

// translate ids to component identifiers
// and add the edges to newgraph

static tpair *translate (matrix M, uint64_t *size)

{ uint r,c,cr,cc;
    uint64_t j,pos,next;
    tpair *newgraph = (tpair*)myalloc(M->elems*sizeof(tpair));
    pos = 0;
    for (r=0;r<M->nrows;r++)
    { cr = info[r].compid;
        next = access64(M->rowpos,r+1,bitsM);
        for (j=access64(M->rowpos,r,bitsM);j<next;j++)
        { c = hashSearch(hash,M->colsbyrow[j]);
            if (c == noval) cc = ncomp+M->colsbyrow[j]; // sink, no row id
            else cc = info[c].compid;
            // heuristic to avoid too many repeated edges
            if (!pos ||
                (newgraph[pos-1].r != cr)||(newgraph[pos-1].c != cc))
            { newgraph[pos].r = cr; newgraph[pos].c = cc;
                pos++;
            }
        }
    }
    *size = pos;
    return newgraph;
}

static matrix strongly (matrix M)

{ uint i,r,c,n,side;
    struct s_matrix C;
    tpair *newgraph;
    uint *edges; // edge pairs
    uint64_t pos,nedges;
    matrix A;

    n = M->nrows; // a node not here has its own connected component
    if (n == 0) return matEmpty(M->height,M->width);

    // creates hash for row values
    hash = hashCreate (M->rowids,M->nrows);

    // initialize the process
    tindex = 1;
    stack = (uint*)myalloc(n*sizeof(uint));
    info = (tarjan*)myalloc(n*sizeof(tarjan));
    compnodes = (uint*)myalloc(n*sizeof(uint));
    comptr = (uint*)myalloc((n+1)*sizeof(uint));
    ncomp = 0; pcomp = 0;
    for (i=0;i<n;i++) info[i].index = 0; // undefined

    // find components
    bitsM = numbits(M->elems);
    for (i=0; i<n ;i++)
        if (info[i].index == 0) strongconnect(M,i);
    myfree(stack);

    comptr[ncomp] = pcomp;
    comptr = (uint*)myrealloc(comptr,(ncomp+1)*sizeof(uint));
    compnodes = (uint*)myrealloc(compnodes,pcomp*sizeof(uint));

    // translate ids to component identifiers
    // and add the edges to newgraph
    newgraph = translate(M,&pos);
    myfree(info);
    hashDestroy(hash);

    // remove possibly many repeated edges and create reduced matrix
    C.colsbyrow = (uint*)myalloc(pos*sizeof(uint));
    C.rowpos = (uint64_t*)myalloc((ncomp+1)*sizeof(uint64_t));
    createReduced (newgraph,pos,C.colsbyrow,C.rowpos);


    C.nrows = ncomp;
    C.elems = C.rowpos[C.nrows];
    myfree(newgraph);
    C.colsbyrow = (uint*)myrealloc(C.colsbyrow,C.elems*sizeof(uint));

    // propagate descendants in topological order
    propagate (C);
    myfree(C.rowpos); myfree(C.colsbyrow);

    // expand to original nodes
    edges = expand (ncomp,&nedges);
    myfree(compnodes);
    myfree(comptr);

    // create matrix from it
    side = mmax(M->width,M->height);
    A = matCreate(side,side,nedges,edges);
    myfree(edges);
    return A;
}

matrix matClos (matrix A, uint pos)

{ matrix M,S,Id;

    S = strongly(A);
    if (!pos)
    { Id = matId(mmax(A->width,A->height));
        M = S; S = matSum(S,Id); matDestroy(M);
        matDestroy(Id);
    }
    return S;
}

// multiplies M by vector V. allocates and returns the result M x V
// use matrix transposition to do V^T x M = M^T x V (transposed vector)
// V is assumed to be of size A->width, the output is of size A->height

double *matVectorMult (matrix A, double *vec)

{ double *W; // W = MxV
    uint i,bits;
    uint *cols,*top;
    double *rowW;

    if (A->height == 0) return NULL;
    W = (double*)myalloc(A->height * sizeof(double));
    for (i=0;i<A->height;i++) W[i] = 0;
    if (!A->elems) return W;
    bits = numbits(A->elems);
    cols = A->colsbyrow;
    for (i=0;i<A->nrows;i++)
    { rowW = W + A->rowids[i];
        top = A->colsbyrow + access64(A->rowpos,i+1,bits);
        while (cols<top) *rowW += vec[*cols++];
    }
    return W;
}

