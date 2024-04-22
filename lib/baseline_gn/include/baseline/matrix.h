
#ifndef INCLUDEDmatrix2
#define INCLUDEDmatrix2

// 32-bit version

#include "basics.h"

typedef struct s_matrix {
    uint elems; // how many filled cells
    uint width,height; // matrix dimensions
    uint nrows; // # nonempty rows
    uint *rowids; // ids of nonempty rows
    uint64_t *rowpos; // starting positions of rows in colsbyrow
    uint *colsbyrow; // cols sorted by row:col
    uint ncols; // etc
    uint *colids;
    uint64_t *colpos;
    uint *rowsbycol;
} *matrix;

// creates matrix of width x height with n cells (2n ints row,col)
// reorders cells array
matrix matCreate (uint height, uint width, uint64_t n, uint *cells);

// destroys matrix
void matDestroy (matrix M);

// creates an empty matrix
matrix matEmpty (uint height, uint width);

// creates an identity matrix
matrix matId (uint side);

// creates a new copy of M, with its own data
matrix matCopy (matrix M);

// transpose a matrix, creating a non-allocated copy that shares the
// data. You need not (and should not) matDestroy this copy
// use *M = matTranspose(M) to actually transpose M
struct s_matrix matTranspose (matrix M);

// writes M to file, which must be opened for writing
void matSave (matrix M, FILE *file);

// loads matrix from file, which must be opened for reading
matrix matLoad (FILE *file);

// space of the matrix, in 32-bit words
uint64_t matSpace (matrix M);

// dimensions of M, returns #elems and writes the others if not null
uint64_t matDims (matrix M, uint *width, uint *height);

// accesses a cell
uint matAccess (matrix M, uint row, uint col);

// recovers all the cells in [r1..r2] x [c1..c2]
// writes 2n integers in buffer, which must be preallocated
// to 2*elems of uint (worst case)
// returns number of elements. just counts if buffer is NULL
uint64_t matCollect (matrix M, uint r1, uint r2,
                     uint c1, uint c2, uint *buffer);

#define fullSide (~(uint)0)

// (boolean) sum of two matrices, assumed to be of the same side
matrix matSum (matrix A, matrix B);

// version with one row or one column, or both
matrix matSum1 (uint row, matrix A, matrix B, uint col);

// (boolean) product of two matrices, assumed to be of the same side
// only rowA of A and colB of B are considered if not fullSide

matrix matMult (matrix A, matrix B);

// version with one row or one column, or both
matrix matMult1 (uint row, matrix A, matrix B, uint col);

// transitive closure of a matrix, pos says if it's + rather than *
matrix matClos (matrix A, uint pos);

// versions to choose one row or one column, or both
matrix matClos1 (uint row, matrix A, uint pos, uint col);

// computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
matrix matMultClos1 (uint row, matrix A, matrix B, uint pos, uint col);

// computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
matrix matClosMult1 (uint row, matrix A, uint pos, matrix B, uint col);

// multiplies M by vector V. allocates and returns the result M x V
// use matrix transposition to do V^T x M = M^T x V (transposed vector)
// V is assumed to be of size A->width, the output is of size A->height
double *matVectorMult (matrix A, double *vec);

#endif
