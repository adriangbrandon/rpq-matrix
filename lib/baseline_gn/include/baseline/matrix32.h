
#ifndef INCLUDEDmatrix32
#define INCLUDEDmatrix32

	// 32-bit version

#include "basics.h"

typedef struct s_matrix32 {
   uint elems; // how many filled cells
   uint width,height; // matrix dimensions
   uint nrows; // # nonempty rows
   uint *rowids; // ids of nonempty rows
   uint *rowpos; // starting positions of rows in colsbyrow
   uint *colsbyrow; // cols sorted by row:col
   uint ncols; // etc
   uint *colids; 
   uint *colpos; 
   uint *rowsbycol;
   } *matrix32;

	// creates matrix of width x height with n cells (2n ints row,col) 
	// reorders cells array
matrix32 matCreate32 (uint height, uint width, uint n, uint *cells);

	// destroys matrix
void matDestroy32 (matrix32 M);
   
        // creates an empty matrix
matrix32 matEmpty32 (uint height, uint width);

        // creates an identity matrix
matrix32 matId32 (uint side);

        // creates a new copy of M, with its own data
matrix32 matCopy32 (matrix32 M);

        // transpose a matrix, creating a non-allocated copy that shares the
        // data. You need not (and should not) matDestroy this copy
        // use *M = matTranspose(M) to actually transpose M
struct s_matrix32 matTranspose32 (matrix32 M);

        // writes M to file, which must be opened for writing
void matSave32 (matrix32 M, FILE *file);

        // loads matrix from file, which must be opened for reading
matrix32 matLoad32 (FILE *file);

        // space of the matrix, in 32-bit words
uint64_t matSpace32 (matrix32 M);

	// dimensions of M, returns #elems and writes the others if not null
uint matDims32 (matrix32 M, uint *width, uint *height);

	// accesses a cell
uint matAccess32 (matrix32 M, uint row, uint col);

	// recovers all the cells in [r1..r2] x [c1..c2]
	// writes 2n integers in buffer, which must be preallocated
	// to 2*elems of uint (worst case)
	// returns number of elements. just counts if buffer is NULL
uint matCollect32 (matrix32 M, uint r1, uint r2,
		 uint c1, uint c2, uint *buffer);

#define fullSide (~(uint)0)

        // (boolean) sum of two matrices, assumed to be of the same side
matrix32 matSum32 (matrix32 A, matrix32 B);

        // version with one row or one column, or both
matrix32 matSum132 (uint row, matrix32 A, matrix32 B, uint col);

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide

matrix32 matMult32 (matrix32 A, matrix32 B);

        // version with one row or one column, or both
matrix32 matMult132 (uint row, matrix32 A, matrix32 B, uint col);

        // transitive closure of a matrix, pos says if it's + rather than *
matrix32 matClos32 (matrix32 A, uint pos);

        // versions to choose one row or one column, or both
matrix32 matClos132 (uint row, matrix32 A, uint pos, uint col);

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
matrix32 matMultClos132 (uint row, matrix32 A, matrix32 B, uint pos, uint col);

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
matrix32 matClosMult132 (uint row, matrix32 A, uint pos, matrix32 B, uint col);

matrix32 matAnd32(matrix32 A, matrix32 B);

matrix32 matAnd132(uint row, matrix32 A, matrix32 B, uint col);

#endif
