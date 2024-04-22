
// 32-bit version

#include <stdlib.h>
#include <string.h>

#include "baseline/matrix32.h"
#include "baseline/heap.h"
#include "baseline/hash.h"
#include "baseline/utilstime.h"
#define w32 32

// creates matrix32 of width x height with n cells (2n ints row,col)
// reorders cells array

// builds row or col

static void sort32(uint* r, int lo, int up )

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
        if ( i-lo < up-i ) { sort32(r,lo,i-1);  lo = i+1; }
        else    { sort32(r,i+1,up);  up = i-1; }
    }
}

static uint buildDim32 (uint size, uint nelems, uint *cells, uint fst, uint snd,
                      uint **ids, uint **pos, uint *other)

{ uint i,d,t;
    uint *dim;

    dim = (uint*)myalloc((1+size)*sizeof(uint));
    for (d=0;d<=size;d++) dim[d] = 0;
    for (i=0;i<nelems;i++) dim[1+cells[2*i+fst]]++;
    t = 0;
    for (d=1;d<=size;d++)
        if (dim[d]) t++;
    *ids = (uint*)myalloc((t+1)*sizeof(uint));
    *pos = (uint*)myalloc((t+1)*sizeof(uint));
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
    { sort32(other,(*pos)[d],(*pos)[d+1]-1);
    }
    return t;
}

matrix32 matCreate32 (uint height, uint width, uint n, uint *cells)

{ matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->height = height;
    M->width = width;
    M->elems = n;
    if (n == 0)
    { M->rowids = M->rowpos = M->colsbyrow = NULL;
        M->colids = M->colpos = M->rowsbycol = NULL;
        M->nrows = M->ncols = 0;
    }
    else
    { M->colsbyrow = (uint*)myalloc(n*sizeof(uint));
        M->rowsbycol = (uint*)myalloc(n*sizeof(uint));
        M->nrows = buildDim32(height,n,cells,0,1,&M->rowids,&M->rowpos,
                            M->colsbyrow);
        M->ncols = buildDim32(width,n,cells,1,0,&M->colids,&M->colpos,
                            M->rowsbycol);
    }
    return M;
}


// destroys matrix32
void matDestroy32 (matrix32 M)

{ if (M == NULL) return;
    myfree (M->rowids); myfree (M->rowpos); myfree (M->colsbyrow);
    myfree (M->colids); myfree (M->colpos); myfree (M->rowsbycol);
    myfree(M);
}

// creates an empty matrix32
matrix32 matEmpty32 (uint height, uint width)

{ matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->height = height;
    M->width = width;
    M->elems = 0;
    M->rowids = M->rowpos = M->colsbyrow = NULL;
    M->colids = M->colpos = M->rowsbycol = NULL;
    M->nrows = M->ncols = 0;
    return M;
}

// creates an identity matrix32
matrix32 matId32 (uint side)

{ matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    uint i;
    M->height = M->width = side;
    M->elems = M->nrows = M->ncols = side;
    M->rowids = (uint*)myalloc((side+1)*sizeof(uint));
    M->rowpos = (uint*)myalloc((side+1)*sizeof(uint));
    M->colids = (uint*)myalloc((side+1)*sizeof(uint));
    M->colpos = (uint*)myalloc((side+1)*sizeof(uint));
    M->colsbyrow = (uint*)myalloc(side*sizeof(uint));
    M->rowsbycol = (uint*)myalloc(side*sizeof(uint));
    for (i=0;i<side;i++)
    { M->rowids[i] = M->rowpos[i] = M->colsbyrow[i] = i;
        M->colids[i] = M->colpos[i] = M->rowsbycol[i] = i;
    }
    M->rowids[side] = M->colids[side] = side;
    M->rowpos[side] = M->colpos[side] = side;
    return M;
}

// creates a new copy of M, with its own data
matrix32 matCopy32 (matrix32 M)

{ matrix32 N = (matrix32)myalloc(sizeof(struct s_matrix32));
    *N = *M;
    if (N->elems > 0)
    { N->rowids = (uint*)myalloc((N->nrows+1)*sizeof(uint));
        memcpy (N->rowids,M->rowids,(N->nrows+1)*sizeof(uint));
        N->rowpos = (uint*)myalloc((N->nrows+1)*sizeof(uint));
        memcpy (N->rowpos,M->rowpos,(N->nrows+1)*sizeof(uint));
        N->colsbyrow = (uint*)myalloc(N->elems*sizeof(uint));
        memcpy(N->colsbyrow,M->colsbyrow,N->elems*sizeof(uint));
        N->colids = (uint*)myalloc((N->ncols+1)*sizeof(uint));
        memcpy (N->colids,M->colids,(N->ncols+1)*sizeof(uint));
        N->colpos = (uint*)myalloc((N->ncols+1)*sizeof(uint));
        memcpy (N->colpos,M->colpos,(N->ncols+1)*sizeof(uint));
        N->rowsbycol = (uint*)myalloc(N->elems*sizeof(uint));
        memcpy(N->rowsbycol,M->rowsbycol,N->elems*sizeof(uint));
    }
    return N;
}

// transpose a matrix32, creating a non-allocated copy that shares the
// data. You need not (and should not) matDestroy this copy
// use *M = matTranspose(M) to actually transpose M


struct s_matrix32 matTranspose32 (matrix32 M)

{ uint aux;
    uint *paux;
    struct s_matrix32 T;
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
void matSave32 (matrix32 M, FILE *file)

{ fwrite(&M->elems,sizeof(uint),1,file);
    fwrite(&M->width,sizeof(uint),1,file);
    fwrite(&M->height,sizeof(uint),1,file);
    if (M->elems == 0) return;
    fwrite(&M->nrows,sizeof(uint),1,file);
    fwrite(&M->ncols,sizeof(uint),1,file);

    fwrite(M->rowids,sizeof(uint),M->nrows+1,file);
    fwrite(M->rowpos,sizeof(uint),M->nrows+1,file);
    fwrite(M->colsbyrow,sizeof(uint),M->elems,file);

    fwrite(M->colids,sizeof(uint),M->ncols+1,file);
    fwrite(M->colpos,sizeof(uint),M->ncols+1,file);
    fwrite(M->rowsbycol,sizeof(uint),M->elems,file);
}

// loads matrix32 from file, which must be opened for reading
matrix32 matLoad32 (FILE *file)

{ matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    fread(&M->elems,sizeof(uint),1,file);
    fread(&M->width,sizeof(uint),1,file);
    fread(&M->height,sizeof(uint),1,file);
    if (M->elems == 0)
    { M->rowids = M->rowpos = M->colsbyrow = NULL;
        M->colids = M->colpos = M->rowsbycol = NULL;
        M->nrows = M->ncols = 0;
    }
    else
    { fread(&M->nrows,sizeof(uint),1,file);
        fread(&M->ncols,sizeof(uint),1,file);

        M->rowids = (uint*)myalloc((M->nrows+1)*sizeof(uint));
        fread(M->rowids,sizeof(uint),(M->nrows+1),file);
        M->rowpos = (uint*)myalloc((M->nrows+1)*sizeof(uint));
        fread(M->rowpos,sizeof(uint),(M->nrows+1),file);
        M->colsbyrow = (uint*)myalloc(M->elems*sizeof(uint));
        fread(M->colsbyrow,sizeof(uint),M->elems,file);

        M->colids = (uint*)myalloc((M->ncols+1)*sizeof(uint));
        fread(M->colids,sizeof(uint),(M->ncols+1),file);
        M->colpos = (uint*)myalloc((M->ncols+1)*sizeof(uint));
        fread(M->colpos,sizeof(uint),(M->ncols+1),file);
        M->rowsbycol = (uint*)myalloc(M->elems*sizeof(uint));
        fread(M->rowsbycol,sizeof(uint),M->elems,file);
    }
    return M;
}

// space of the matrix32, in 32-bit words

uint64_t matSpace32 (matrix32 M)

{ return (sizeof(struct s_matrix32)/sizeof(uint) +
          2 * (M->nrows+1) + 2 * (M->ncols+1) + 2 * (uint64_t)M->elems)
         * (sizeof(uint) / 4);
}

// dimensions of M, returns #elems and writes the others if not null
uint matDims32 (matrix32 M, uint *width, uint *height)

{ if (width) *width = M->width;
    if (height) *height = M->height;
    return M->elems;
}

// accesses a cell

// finds first >= key in r[low..high], may return high+1
static int search32(uint key, int low, int high, uint *r)

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

uint matAccess32 (matrix32 M, uint row, uint col)

{  // search by row, arbitrarily
    int rpos,rcol;
    if (M->elems == 0) return 0;
    rpos = search32(row,0,M->nrows-1,M->rowids);
    if ((rpos == M->nrows) || (M->rowids[rpos] != row)) return 0;
    rcol = search32(col,M->rowpos[rpos],M->rowpos[rpos+1]-1,M->colsbyrow);
    if ((rcol == M->rowpos[rpos+1]) || (M->colsbyrow[rcol] != col)) return 0;
    return 1;
}

// recovers all the cells in [r1..r2] x [c1..c2]
// writes 2n integers in buffer, which must be preallocated
// to 2*elems of uint (worst case)
// returns number of elements. just counts if buffer is NULL

uint matCollect32 (matrix32 M, uint r1, uint r2,
                 uint c1, uint c2, uint *buffer)

{ int rpos,rcol;
    uint nrows,rowid,erow;
    uint t = 0;
    if (M->elems == 0) return 0;
    nrows = M->nrows;
    rpos = search32(r1,0,nrows-1,M->rowids);
    while ((rpos < nrows) && ((rowid = M->rowids[rpos]) <= r2))
    { erow = M->rowpos[rpos+1];
        rcol = search32(c1,M->rowpos[rpos],erow-1,M->colsbyrow);
        if (buffer)
        { while ((rcol < erow) && (M->colsbyrow[rcol] <= c2))
            { buffer[2*t] = rowid;
                buffer[2*t+1] = M->colsbyrow[rcol];
                t++;
                rcol++;
            }
        }
        else if (rcol < erow)
        { t -= rcol;
            rcol = search32(c2+1,rcol,erow-1,M->colsbyrow);
            t += rcol;
        }
        rpos++;
    }
    return t;
}

// (boolean) sum of two matrices, assumed to be of the same side

static uint merge32 (uint *dataA, uint fa, uint *dataB, uint fb, uint *dataM)

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

static uint matMerge32 (uint na, uint nb, uint *posA, uint *posB,
                      uint *posM, uint *idsA, uint *idsB, uint *idsM,
                      uint *dataA, uint *dataB, uint *dataM)

{ uint c,p,pa,pb,dif;
    c = p = pa = pb = 0;
    while ((pa < na) || (pb < nb))
    { posM[p] = c;
        //if (pa < na && (pb == nb || idsA[pa] < idsB[pb]))
        if (idsA[pa] < idsB[pb])
        { idsM[p] = idsA[pa];
            dif = posA[pa+1]-posA[pa];
            memcpy(dataM+c,dataA+posA[pa],dif*sizeof(uint));
            pa++; c += dif;
        }
        //else if (pb < nb && (pa == na || idsB[pb] < idsA[pa]))
        else if (idsB[pb] < idsA[pa])
        { idsM[p] = idsB[pb];
            dif = posB[pb+1]-posB[pb];
            memcpy(dataM+c,dataB+posB[pb],dif*sizeof(uint));
            pb++; c += dif;
        }
        else
        { idsM[p] = idsA[pa];
            c += merge32(dataA+posA[pa],posA[pa+1]-posA[pa],
                       dataB+posB[pb],posB[pb+1]-posB[pb],dataM+c);
            pa++; pb++;
        }
        p++;
    }
    posM[p] = c;
    return p;
}

matrix32 matSum32 (matrix32 A, matrix32 B)

{ matrix32 M;
    uint p,pa,pb,na,nb,c;

    if ((A->elems == 0) || (B->elems == 0))
    { if (A->elems == 0) M = matCopy32(B); else M = matCopy32(A);
        M->width = mmax(A->width,B->width);
        M->height = mmax(A->height,B->height);
        return M;
    }
    M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->width = mmax(A->width,B->width);
    M->height = mmax(A->height,B->height);

    M->rowids = (uint*)myalloc((A->nrows+B->nrows+1)*sizeof(uint));
    M->rowpos = (uint*)myalloc((A->nrows+B->nrows+1)*sizeof(uint));
    M->colsbyrow = (uint*)myalloc((A->elems+B->elems)*sizeof(uint));
    M->nrows = matMerge32(A->nrows,B->nrows,A->rowpos,B->rowpos,M->rowpos,
                         A->rowids,B->rowids,M->rowids,
                         A->colsbyrow,B->colsbyrow,M->colsbyrow);
    M->rowids = (uint*)myrealloc(M->rowids,(M->nrows+1)*sizeof(uint));
    M->rowpos = (uint*)myrealloc(M->rowpos,(M->nrows+1)*sizeof(uint));
    M->elems = M->rowpos[M->nrows];
    M->rowids[M->nrows] = M->height;
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,M->elems*sizeof(uint));

    M->colids = (uint*)myalloc((A->ncols+B->ncols+1)*sizeof(uint));
    M->colpos = (uint*)myalloc((A->ncols+B->ncols+1)*sizeof(uint));
    M->rowsbycol = (uint*)myalloc(M->elems*sizeof(uint));
    M->ncols = matMerge32 (A->ncols,B->ncols,A->colpos,B->colpos,M->colpos,
                         A->colids,B->colids,M->colids,
                         A->rowsbycol,B->rowsbycol,M->rowsbycol);
    M->colids = (uint*)myrealloc(M->colids,(M->ncols+1)*sizeof(uint));
    M->colpos = (uint*)myrealloc(M->colpos,(M->ncols+1)*sizeof(uint));
    M->colids[M->ncols] = M->width;

    return M;
}

// version with one row or one column, or both

static matrix32 matOne32 (uint height, uint width, uint row, uint col)

{ matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->elems = M->nrows = M->ncols = 1;
    M->height = height; M->width = width;

    M->rowids = (uint*)myalloc(2*sizeof(uint));
    M->rowids[0] = row; M->rowids[1] = height;
    M->rowpos = (uint*)myalloc(2*sizeof(uint));
    M->rowpos[0] = 0; M->rowpos[1] = 1;
    M->colsbyrow = (uint*)myalloc(1*sizeof(uint));
    M->colsbyrow[0] = col;

    M->colids = (uint*)myalloc(2*sizeof(uint));
    M->colids[0] = col; M->colids[1] = width;
    M->colpos = (uint*)myalloc(2*sizeof(uint));
    M->colpos[0] = 0; M->colpos[1] = 1;
    M->rowsbycol = (uint*)myalloc(1*sizeof(uint));
    M->rowsbycol[0] = row;

    return M;
}

static void copy132 (uint **colids, uint **colpos, uint **rowids, uint **rowpos,
                   uint **rowsbycol, uint **colsbyrow, uint *data,
                   uint col, uint width, uint height, uint nelems)

{ uint i;
    *colids = (uint*)myalloc(2*sizeof(uint));
    *colpos = (uint*)myalloc(2*sizeof(uint));
    (*colids)[0] = col; (*colids)[1] = width;
    (*colpos)[0] = 0; (*colpos)[1] = nelems;
    *rowsbycol = (uint*)myalloc(nelems*sizeof(uint));
    memcpy(*rowsbycol,data,nelems*sizeof(uint));
    *rowids = (uint*)myalloc((nelems+1)*sizeof(uint));
    *rowpos = (uint*)myalloc((nelems+1)*sizeof(uint));
    memcpy(*rowids,data,nelems*sizeof(uint));
    (*rowids)[nelems] = height;
    *colsbyrow = (uint*)myalloc(nelems*sizeof(uint));
    for (i=0;i<=nelems;i++) (*rowpos)[i] = i;
    for (i=0;i<nelems;i++) (*colsbyrow)[i] = col;
}

matrix32 matSum132 (uint row, matrix32 A, matrix32 B, uint col)

{ uint posA,posB,i;
    uint height = mmax(A->height,B->height);
    uint width = mmax(A->width,B->width);
    matrix32 M;
    struct s_matrix32 At,Bt;

    if (row == fullSide)
        if (col == fullSide)
            return matSum32(A,B);
        else
        { posA = search32(col,0,A->ncols-1,A->colids);
            posB = search32(col,0,B->ncols-1,B->colids);
            if (((posA >= A->ncols) || (A->colids[posA] != col)) &&
                ((posB >= B->ncols) || (B->colids[posB] != col)))
                return matEmpty32 (height,width);
            M = (matrix32)myalloc(sizeof(struct s_matrix32));
            M->height = height; M->width = width;
            M->ncols = 1;
            if ((posA >= A->ncols) || (A->colids[posA] != col))
            { M->elems = B->colpos[posB+1]-B->colpos[posB];
                copy132 (&M->colids,&M->colpos,&M->rowids,&M->rowpos,
                       &M->rowsbycol,&M->colsbyrow,
                       B->rowsbycol+B->colpos[posB],
                       col,width,height,M->elems);
            }
            else if ((posB >= B->ncols) || (B->colids[posB] != col))
            { M->elems = A->colpos[posA+1]-A->colpos[posA];
                copy132 (&M->colids,&M->colpos,&M->rowids,&M->rowpos,
                       &M->rowsbycol,&M->colsbyrow,
                       A->rowsbycol+A->colpos[posA],
                       col,width,height,M->elems);
            }
            else
            { uint na = A->colpos[posA+1]-A->colpos[posA];
                uint nb = B->colpos[posB+1]-B->colpos[posB];
                uint *data = (uint*)myalloc((na+nb)*sizeof(uint));
                M->elems = merge32(A->rowsbycol+A->colpos[posA],na,
                                 B->rowsbycol+B->colpos[posB],nb,data);
                copy132 (&M->colids,&M->colpos,&M->rowids,&M->rowpos,
                       &M->rowsbycol,&M->colsbyrow,data,
                       col,width,height,M->elems);
                myfree(data);
            }
            M->nrows = M->elems;
            return M;
        }
    else if (col == fullSide) // sorry for this, I couldn't resist
    { At = matTranspose32(A); Bt = matTranspose32(B);
        M = matSum132(col,&At,&Bt,row);
        *M = matTranspose32(M);
        return M;
    }
    else // both row and col restrictions
    { if (matAccess32(A,row,col) && matAccess32(B,row,col))
            return matOne32 (height,width,row,col);
        else return matEmpty32 (height,width);
    }
}

// (boolean) product of two matrices, assumed to be of the same side
// only row of A and col of B are considered if not fullSide

// assumes na << nb and data is preallocated (presumably to na elems)

static uint inter120_32 (uint *dataA, uint na, uint *dataB, uint nb)

{ int i,j;
    // uint p = 0;
    j = 0;
    for (i=0;i<na;i++)
    { uint d = dataA[i];
        uint s = 1;
        j--;
        while ((j+s < nb) && (dataB[j+s] < d)) s <<= 1;
        j = search32(d,j+1+(s>>1),mmin(j+s,nb-1),dataB);
        if ((j < nb) && (dataB[j] == d)) return 1; // data[p++] = d;
    }
    return 0; // p;
}

static uint inter0_32 (uint *dataA, uint na, uint *dataB, uint nb)

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

static uint inters0_32 (uint *dataA, uint na, uint *dataB, uint nb)

{ if (nb > 30*na) return inter120_32(dataA,na,dataB,nb);
    if (na > 30*nb) return inter120_32(dataB,nb,dataA,na);
    return inter0_32(dataA,na,dataB,nb);
}

typedef struct
{ uint rowB,colA;
    uint ncols,nrows;
    uint *cbyr,*rbyc;
    uint pcbyr;
} ttask;

// assumes task is preallocated, to at least min(na,nb) elements

static uint inter12_32 (uint *dataA, uint na, uint *dataB, uint nb,
                     ttask *task, uint aisr)

{
    int i,j;
    uint p = 0;
    j = 0;
    for (i=0;i<na;i++)
    { uint d = dataA[i];
        uint s = 1;
        j--;
        while ((j+s < nb) && (dataB[j+s] < d)) s <<= 1;
        j = search32(d,j+1+(s>>1),mmin(j+s,nb-1),dataB);
        if ((j < nb) && (dataB[j] == d))
        { if (aisr) { task[p].colA = i; task[p].rowB = j; }
            else { task[p].colA = j; task[p].rowB = i; }
            p++; j++;
        }
    }
    return p;
}

static uint inter_32 (uint *dataA, uint na, uint *dataB, uint nb, ttask *task)

{
    //fprintf(stdout, "Exp\n");
    int i,j;
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


static uint inters_32 (uint *dataA, uint na, uint *dataB, uint nb, ttask *task)

{ if (nb > 30*na) return inter12_32(dataA,na,dataB,nb,task,1);
    if (na > 30*nb) return inter12_32(dataB,nb,dataA,na,task,0);
    return inter_32(dataA,na,dataB,nb,task);
}




matrix32 matMult32 (matrix32 A, matrix32 B)

{

#if TIMEOUT
    if(time_diff() > TIMEOUT) {
        //fprintf(stdout, "Time %llu\n", user_diff());
        return matEmpty32(A->height, B->width);
    }
#endif

    ttask *task;
    uint c,i,j,m,p,pr,pc,prc;
    heap Hr,Hc;
    uint size;
    uint *colc; // column counter to later build the transpose
    uint *stackp,*stackv; // implement colc as an initializable array
    uint stacks;

    // tasks are columns of A that are also rows of B
    task = (ttask*)myalloc(mmin(A->ncols,B->nrows)*sizeof(ttask));
    p = inters_32(A->colids,A->ncols,B->rowids,B->nrows,task);
    if (p == 0) // empty result!
    { myfree(task);
        return matEmpty32(A->height,B->width);
    }
    task = (ttask*)myrealloc(task,p*sizeof(ttask));
    Hr = (heap)myalloc(p*sizeof(struct s_heap)); // to merge rows
    Hc = (heap)myalloc(p*sizeof(struct s_heap)); // to merge columns

    // we first fill M in row-wise form
    matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->rowids = (uint*)myalloc((A->nrows+1)*sizeof(uint));
    M->rowpos = (uint*)myalloc((A->nrows+1)*sizeof(uint));
    size = 1 << numbits(A->nrows);
    M->colsbyrow = (uint*)myalloc(size*sizeof(uint)); // growing array
    pr = prc = 0;
    //colc = (uint*)myalloc(B->width*sizeof(uint));
    colc = (uint*)myalloc(B->width*sizeof(uint));
    stackp = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stackv = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stacks = 0;
    // for (c=0;c<B->width;c++) colc[c] = 0;

    // create the p tasks
    for (i=0;i<p;i++)
    { task[i].ncols = B->rowpos[task[i].rowB+1]-B->rowpos[task[i].rowB];
        task[i].cbyr = B->colsbyrow+B->rowpos[task[i].rowB];
        task[i].nrows = A->colpos[task[i].colA+1]-A->colpos[task[i].colA];
        task[i].rbyc = A->rowsbycol+A->colpos[task[i].colA];
        Hr[i].key = *task[i].rbyc++;
        task[i].nrows--;
        Hr[i].data = i;
    }

    // traverse them by increasing row
    heapify(Hr,p);
    while (p)
    { struct s_heap hr;
        // get all equal minimum hr.key and collect hr.data in Hc
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
	  return matEmpty32(A->height,B->width);
	}
     M->nrows = pr;
     M->rowids[pr] = A->height;
     M->rowpos[pr] = prc;
     M->rowids = (uint*)myrealloc(M->rowids,(pr+1)*sizeof(uint));
     M->rowpos = (uint*)myrealloc(M->rowpos,(pr+1)*sizeof(uint));
     M->colsbyrow = (uint*)myrealloc(M->colsbyrow,prc*sizeof(uint));
     M->elems = prc;

	// now create the transposed representation

     //M->colids = (uint*)myalloc((B->width+1)*sizeof(uint));
     //pc = 0;
     // for (c=0;c<B->width;c++)
     //     { if (colc[c]) M->colids[pc++] = c; }
    M->colids = stackp;
    for (pc=0;pc<stacks;pc++){
        colc[stackp[pc]] = stackv[pc];
    }
    sort32(M->colids,0,pc-1);
     M->colids[pc] = B->width;
     M->ncols = pc;
     M->colids = (uint*)myrealloc(M->colids,(pc+1)*sizeof(uint));
     M->colpos = (uint*)myalloc((pc+1)*sizeof(uint));
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

    M->width = B->width;
    M->height = A->height;
    return M;
}

matrix32 matMult1R32 (uint r, matrix32 A, matrix32 B)

{ ttask *task;
    uint i,j,p,pr,pc,prc;
    heap Hc;
    matrix32 M;

    pr = search32(r,0,A->nrows,A->rowids);
    if (A->rowids[pr] != r) return matEmpty32(A->height,B->width);

    // tasks are columns of A that are also rows of B
    task = (ttask*)myalloc(mmin(A->ncols,B->nrows)*sizeof(ttask));
    p = inters_32(A->colsbyrow+A->rowpos[pr],A->rowpos[pr+1]-A->rowpos[pr],
               B->rowids,B->nrows,task);
    if (p == 0) // empty result!
    { myfree(task);
        return matEmpty32(A->height,B->width);
    }
    task = (ttask*)myrealloc(task,p*sizeof(ttask));
    Hc = (heap)myalloc(p*sizeof(struct s_heap)); // to merge columns

    M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->colsbyrow = (uint*)myalloc(B->ncols*sizeof(uint));
    prc = 0;

    // create the p tasks
    for (i=0;i<p;i++)
    { task[i].ncols = B->rowpos[task[i].rowB+1]-B->rowpos[task[i].rowB];
        task[i].cbyr = B->colsbyrow+B->rowpos[task[i].rowB];
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
        return matEmpty32(A->height,B->width);
    }
    M->nrows = 1;
    M->rowids = (uint*)myalloc(2*sizeof(uint));
    M->rowids[0] = r; M->rowids[1] = A->height;
    M->rowpos = (uint*)myalloc(2*sizeof(uint));
    M->rowpos[0] = 0; M->rowpos[1] = prc;
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,prc*sizeof(uint));

    // create the transposed version
    M->ncols = prc;
    M->colids = (uint*)myalloc((prc+1)*sizeof(uint));
    M->colpos = (uint*)myalloc((prc+1)*sizeof(uint));
    M->rowsbycol = (uint*)myalloc(prc*sizeof(uint));
    memcpy(M->colids,M->colsbyrow,prc*sizeof(uint));
    for (pc=0;pc<prc;pc++)
    { M->colpos[pc] = pc;
        M->rowsbycol[pc] = r;
    }
    M->colpos[pc] = pc;
    M->colids[pc] = B->width;

    M->width = B->width;
    M->height = A->height;
    M->elems = prc;
    return M;
}

// version with one row or one column, or both

matrix32 matMult132 (uint row, matrix32 A, matrix32 B, uint col)

{ matrix32 M;
    struct s_matrix32 At,Bt;

    if ((A->elems == 0) || (B->elems == 0))
        return matEmpty32(A->height,B->width);
    if ((row == fullSide) && (col == fullSide)) return matMult32(A,B);
    if ((row != fullSide) && (col != fullSide))
    { if (inters0_32(A->colsbyrow+A->rowpos[row],
                  A->rowpos[row+1]-A->rowpos[row],
                  B->rowsbycol+B->colpos[col],
                  B->colpos[col+1]-B->colpos[col]))
            return matOne32(A->height,B->width,row,col);
        return matEmpty32(A->height,B->width);
    }
    if (row == fullSide)
    { At = matTranspose32(A); Bt = matTranspose32(B);
        M = matMult1R32(col,&Bt,&At);
        *M = matTranspose32(M);
        return M;
    }
    return matMult1R32(row,A,B);
}

// transitive closure of a matrix32, pos says if it's + rather than *

matrix32 matClos032 (matrix32 A, uint pos)

{ matrix32 M,P,S,Id;
    uint elems;

    elems = A->elems;
    P = matMult32 (A,A);
    S = matSum32 (A,P);
    while (S->elems != elems)
    { elems = S->elems;
        matDestroy32(P);
        P = matMult32(S,S);
        M = matSum32(S,P);
        matDestroy32(S);
        S = M;
    }
    matDestroy32(P);
    if (!pos)
    { Id = matId32(mmax(A->width,A->height));
        M = S; S = matSum32(S,Id); matDestroy32(M);
        matDestroy32(Id);
    }
    return S;
}

// versions to chose row row or column col at the end
// here we restrict first to the row/column even if we lose the
// exponential increase in the paths

// if coltest not null, checks that the cell *coltest exists and
// returns there 1 or 0 as soon as it can
// if not null, I implies starting with I, I x A* or I x A+

static matrix32 matClosRow32 (uint row, matrix32 I, matrix32 A, uint pos, uint *coltest)

{ matrix32 M,P,S,E;
    uint elems;
    uint dim = mmax(A->height,A->width);

    if (I == NULL)
    { if (pos) E = matEmpty32(dim,dim);
        else E = matOne32(dim,dim,row,row);
        S = matSum132 (row,A,E,fullSide);
        matDestroy32(E);
    }
    else
    { if (pos) S = matMult132(row,I,A,fullSide);
        else { E = matEmpty32(dim,dim);
            S = matSum132(row,I,E,fullSide);
            matDestroy32(E);
        }
    }
    elems = S->elems;
    if (coltest && matAccess32(S,row,*coltest))
    { matDestroy32(S); *coltest = 1; return NULL; }
    P = matMult32 (S,A);
    M = S; S = matSum32 (P,S); matDestroy32(M);
    while (S->elems != elems)
    { if (coltest && matAccess32(S,row,*coltest))
        { matDestroy32(S); matDestroy32(P); *coltest = 1; return NULL; }
        elems = S->elems;
        M = P; P = matMult32(P,A); matDestroy32(M);
        M = S; S = matSum32(S,P); matDestroy32(M);
    }
    matDestroy32(P);
    if (coltest) { matDestroy32(S); *coltest = 0; return NULL; }
    return S;
}

// versions to choose one row or one column, or both

matrix32 matClos132 (uint row, matrix32 A, uint pos, uint col)

{ uint nrow,ncol;
    uint side = mmax(A->width,A->height);
    uint test;
    matrix32 M;
    struct s_matrix32 At;

    if (row == fullSide)
    { if (col == fullSide) return matClos32(A,pos);
        else
        { At = matTranspose32(A);
            M = matClosRow32(col,NULL,&At,pos,NULL);
            *M = matTranspose32(M);
            return M;
        }
    }
    else
    { if (col == fullSide) return matClosRow32(row,NULL,A,pos,NULL);
    }
    // both row and col
    nrow = matCollect32 (A,0,fullSide,col,col,NULL);
    ncol = matCollect32 (A,row,row,0,fullSide,NULL);
    if (ncol < nrow)
    { test = row;
        At = matTranspose32(A);
        matClosRow32(col,NULL,&At,pos,&test); // does not return matrix32
    }
    else
    { test = col;
        matClosRow32(row,NULL,A,pos,&test); // does not return matrix32
    }
    if (test) return matOne32 (side,side,row,col);
    else return matEmpty32(side,side);
}

// computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)

matrix32 matMultClos132 (uint row, matrix32 A, matrix32 B, uint pos, uint col)

{ uint nrow,ncol;
    uint side = mmax(A->width,A->height);
    uint test;
    matrix32 M,M1;
    if (row == fullSide)
    { if (col == fullSide)  // nothing special
        { M1 = matClos32(B,pos);
            M = matMult32(A,M1);
            matDestroy32(M1);
            return M;
        }
        else // it's much better to start with restricted B
        { M1 = matClos132(fullSide,B,pos,col);
            M = matMult32(A,M1);
            matDestroy32(M1);
            return M;
        }
    }
    else
    { if (col == fullSide) return matClosRow32(row,A,B,pos,NULL);
    }
    // both row and col
    test = col;
    matClosRow32(row,A,B,pos,&test); // returns no matrix32, just test
    if (test) return matOne32 (side,side,row,col);
    return matEmpty32(side,side);
}

// computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)

matrix32 matClosMult132 (uint row, matrix32 A, uint pos, matrix32 B, uint col)

{ matrix32 M;
    struct s_matrix32 At,Bt;
    At = matTranspose32(A);
    Bt = matTranspose32(B);
    M = matMultClos132(col,&Bt,&At,pos,row);
    *M = matTranspose32(M);
    return M;
}


// finds the strongly connected components of the matrix32

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

static void strongconnect (matrix32 M, uint v)

{ uint c,r,u;
    info[v].index  = tindex;
    info[v].lowlink = tindex;
    info[v].compid = onstack;
    tindex++;
    stack[pstack++] = v;

    for (c=M->rowpos[v]; c<M->rowpos[v+1]; c++)
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

static void createReduced (tpair *H, uint n,
                           uint *colsbyrow, uint *rowpos)

{ uint i;
    uint p,b; // writes at colsbyrow[p], heap is H[b..n-1]
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

static inline void siftdown (uint *H, uint n, uint i)

{ uint val = H[n-i];
    uint pos;
    while (2*i <= n)
    { if ((2*i+1 > n) || (H[n-2*i] < H[n-(2*i+1)])) pos = 2*i;
        else pos = 2*i+1;
        if (H[n-pos] < val) { H[n-i] = H[n-pos]; i = pos; }
        else break;
    }
    H[n-i] = val;
}

static uint heapsort (uint *H, uint n)

{ uint i;
    uint p,b; // writes at p, heap is H[b..n-1]
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
static uint *rowpos;

#define resize(array,size,bits,newsize) \
   { if ((newsize) > size) \
      { if ((newsize) >> bits) \
	   { bits = numbits(newsize); \
             array =(uint*)myrealloc(array,(((uint)1)<<bits)*sizeof(uint)); \
	   } \
	size = (newsize); \
      } \
   }

static void propagate (struct s_matrix32 C)

{ uint size,pos,lsize,bits;
    uint v,u,j;

    size = C.elems;
    bits = numbits(size);
    colsbyrow = (uint*)myalloc((((uint)1)<<bits)*sizeof(uint));
    rowpos = (uint*)myalloc((C.nrows+1)*sizeof(uint));
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
        pos = rowpos[v] + heapsort(colsbyrow+rowpos[v],pos-rowpos[v]);
    }
    rowpos[C.nrows] = pos;
    colsbyrow = (uint*)myrealloc(colsbyrow,pos*sizeof(uint));
}

// expand representatives of components to all their nodes

static uint *expand (uint nrows, uint *nedges)

{ uint *edges;
    uint size,pos,usize,vsize;
    uint v,u,i,j,k;
    uint bits;

    size = 2*rowpos[nrows];
    bits = numbits(size);
    edges = myalloc((((uint)1)<<bits)*sizeof(uint));
    pos = 0;
    for (v=0;v<nrows;v++)
    {  // add product of the connected components
        vsize = comptr[v+1]-comptr[v];
        for (k=rowpos[v];k<rowpos[v+1];k++)
        { u = colsbyrow[k];
            if (u < ncomp)
            { usize = comptr[u+1]-comptr[u];
                resize(edges,size,bits,pos+2*vsize*usize);
                for (i=comptr[v];i<comptr[v+1];i++)
                    for (j=comptr[u];j<comptr[u+1];j++)
                    { edges[pos++] = compnodes[i];
                        edges[pos++] = compnodes[j];
                    }
            }
            else // a sink, real id is u-ncomp
            { resize(edges,size,bits,pos+2*vsize);
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

static tpair *translate (matrix32 M, uint *size)

{ uint j,pos,r,c,cr,cc;
    tpair *newgraph = (tpair*)myalloc(M->elems*sizeof(tpair));
    pos = 0;
    for (r=0;r<M->nrows;r++)
    { cr = info[r].compid;
        for (j=M->rowpos[r];j<M->rowpos[r+1];j++)
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

static matrix32 strongly (matrix32 M)

{ uint i,r,c,n,side;
    struct s_matrix32 C;
    tpair *newgraph;
    uint *edges; // edge pairs
    uint pos,nedges;
    matrix32 A;

    n = M->nrows; // a node not here has its own connected component
    if (n == 0) return matEmpty32(M->height,M->width);

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
    for (i=0;i<n;i++)
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

    // remove possibly many repeated edges and create reduced matrix32
    C.colsbyrow = (uint*)myalloc(pos*sizeof(uint));
    C.rowpos = (uint*)myalloc((ncomp+1)*sizeof(uint));
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

    // create matrix32 from it
    side = mmax(M->width,M->height);
    A = matCreate32(side,side,nedges,edges);
    myfree(edges);
    return A;
}

matrix32 matClos32 (matrix32 A, uint pos)

{ matrix32 M,S,Id;

    S = strongly(A);
    if (!pos)
    { Id = matId32(mmax(A->width,A->height));
        M = S; S = matSum32(S,Id); matDestroy32(M);
        matDestroy32(Id);
    }
    return S;
}


matrix32 matAnd32(matrix32 A, matrix32 B){

    uint c,i,j,p,pc,pr;
    uint size;
    uint *colc; // column counter to later build the transpose
    uint *ptrA, *ptrB;
    uint64_t pos,prc;
    uint *stackp,*stackv; // implement colc as an initializable array
    uint stacks;
    matrix32 M = (matrix32)myalloc(sizeof(struct s_matrix32));
    M->rowids = (uint*)myalloc((A->nrows+1)*sizeof(uint));
    M->rowpos = (uint*)myalloc((A->nrows+1)*sizeof(uint));
    size = 1 << numbits(A->nrows);
    M->colsbyrow = (uint*)myalloc(size*sizeof(uint)); // growing array
    i = j = pc = pr = p = prc = 0;
    colc = (uint*)myalloc(B->width*sizeof(uint));
    stackp = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stackv = (uint*)myalloc((B->ncols+1)*sizeof(uint));
    stacks = 0;
    // for (c=0;c<B->width;c++) colc[c] = 0;

    ttask* task_row, *task_col;
    task_row = (ttask*)myalloc(A->nrows * sizeof(ttask));
    pr = inters_32(A->rowids, A->nrows, B->rowids, B->nrows, task_row);
    if (pr == 0) // empty result!
    { myfree(task_row);
        return matEmpty32(A->height,B->width);
    }
    task_row = (ttask*)myrealloc(task_row, pr * sizeof(ttask));
    while(i < pr){

        task_col = (ttask*)myalloc(A->nrows * sizeof(ttask));
        ptrA = A->colsbyrow + A->rowpos[task_row[i].colA];
        ptrB = B->colsbyrow + B->rowpos[task_row[i].rowB];
        pc = inters_32(ptrA, A->rowpos[task_row[i].colA+1]-A->rowpos[task_row[i].colA],
               ptrB, B->rowpos[task_row[i].rowB+1]-B->rowpos[task_row[i].rowB],
               task_col);
        if(pc > 0){
            M->rowids[p] = B->rowids[task_row[i].rowB];
            M->rowpos[p] = prc;
            ++p;
            task_col = (ttask*)myrealloc(task_col, pc * sizeof(ttask));
            j = 0;
            while(j < pc){
                if (prc == size) // realloc colsbyrow
                { size <<= 1;
                    M->colsbyrow = (uint*)myrealloc(M->colsbyrow, size*sizeof(uint));
                }
                uint k = task_col[j].colA;
                M->colsbyrow[prc++] = k;
                if ((colc[k] >= stacks) || (stackp[colc[k]] != k)) // invalid value
                {   stackp[stacks] = k;
                    stackv[stacks] = 1;
                    colc[k] = stacks++;
                }
                else stackv[colc[k]]++;
                ++j;
            }
        }
        myfree(task_col);
        ++i;
    }

    // close this stage
    myfree(task_row);
    if (prc == 0) // empty result!
    { myfree(M->rowids); myfree(M->rowpos); myfree(M->colsbyrow); myfree(M);
        myfree(colc);myfree(stackv); myfree(stackp);
        return matEmpty32(A->height,B->width);
    }
    M->nrows = p;
    M->rowids[p] = A->height;
    M->rowpos[p] = prc;
    M->rowids = (uint*)myrealloc(M->rowids,(p+1)*sizeof(uint));
    M->rowpos = (uint*)myrealloc(M->rowpos,(p+1)*sizeof(uint));
    M->colsbyrow = (uint*)myrealloc(M->colsbyrow,prc*sizeof(uint));
    M->elems = prc;

    // now create the transposed representation

    M->colids = stackp;
    for (pc=0;pc<stacks;pc++){
        colc[stackp[pc]] = stackv[pc];
    }
    sort32(M->colids,0,pc-1);
    M->colids[pc] = B->width;
    M->ncols = pc;
    M->colids = (uint*)myrealloc(M->colids,(pc+1)*sizeof(uint));
    M->colpos = (uint*)myalloc((pc+1)*sizeof(uint));
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

    M->width = B->width;
    M->height = A->height;
    return M;

}

matrix32 matAnd132(uint row, matrix32 A, matrix32 B, uint col){


    struct s_matrix32 At,Bt;
    matrix32 M;
    if (col == fullSide)
        if (row == fullSide)
            return matAnd32(A,B);
        else {
            uint ra, rb;
            ra = search32(row, 0, A->nrows - 1, A->rowids);
            if ((ra == A->nrows) || (A->rowids[ra] != row)) return matEmpty32 (A->height,B->width);
            rb = search32(row, 0, B->nrows - 1, B->rowids);
            if ((rb == B->nrows) || (B->rowids[rb] != row)) return matEmpty32 (A->height,B->width);

            ttask *task_col;
            uint *ptrA, *ptrB;
            uint pc, j, prc, p, size, c;
            task_col = (ttask *) myalloc(A->nrows * sizeof(ttask));
            ptrA = A->colsbyrow + A->rowpos[ra];
            ptrB = B->colsbyrow + B->rowpos[rb];
            j = pc = p = prc = 0;
            pc = inters_32(ptrA, A->rowpos[ra + 1] - A->rowpos[ra],
                           ptrB, B->rowpos[rb + 1] - B->rowpos[rb],
                           task_col);

            M = (matrix32) myalloc(sizeof(struct s_matrix32));
            if (pc > 0) {
                M->rowids = (uint *) myalloc((2) * sizeof(uint));
                M->rowpos = (uint *) myalloc((2) * sizeof(uint));
                size = 1 << numbits(A->nrows);
                M->colsbyrow = (uint *) myalloc(size * sizeof(uint)); // growing array
                M->rowids[p] = row;
                M->rowpos[p] = prc;
                ++p;
                task_col = (ttask *) myrealloc(task_col, pc * sizeof(ttask));
                j = 0;
                while (j < pc) {
                    if (prc == size) // realloc colsbyrow
                    {
                        size <<= 1;
                        M->colsbyrow = (uint *) myrealloc(M->colsbyrow, size * sizeof(uint));
                    }
                    M->colsbyrow[prc++] = ptrA[task_col[j].colA];
                    ++j;
                }

            }
            myfree(task_col);

            // close this stage
            if (prc == 0) // empty result!
            {
                myfree(M);
                return matEmpty32(A->height, B->width);
            }
            M->nrows = 1;
            M->rowids[p] = M->height;
            M->rowpos[p] = prc;
            M->colsbyrow = (uint *) myrealloc(M->colsbyrow, prc * sizeof(uint));
            M->elems = prc;

            // now create the transposed representation
            M->colids = (uint *) myalloc((prc + 1) * sizeof(uint));
            M->ncols = prc;
            M->colpos = (uint *) myalloc((prc + 1) * sizeof(uint));
            M->rowsbycol = (uint *) myalloc(prc * sizeof(uint));
            for (c = 0; c < prc; c++) {
                M->colpos[c] = c;
                M->colids[c] = M->colsbyrow[c];
                M->rowsbycol[c] = row;
            }
            M->colpos[prc] = prc;
            M->colids[prc] = B->width;
            M->width = B->width;
            M->height = A->height;
            return M;
        }
    else if (row == fullSide) // sorry for this, I couldn't resist
    {
        At = matTranspose32(A); Bt = matTranspose32(B);
        M = matAnd132(col,&At,&Bt,row);
        *M = matTranspose32(M);
        return M;
    }
    else // both row and col restrictions
    { if (matAccess32(A,row,col) && matAccess32(B,row,col))
            return matOne32 (A->height,B->width,row,col);
        else return matEmpty32 (A->height,B->width);
    }

}



