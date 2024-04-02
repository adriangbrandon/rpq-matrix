
#ifndef INCLUDEDmatrix
#define INCLUDEDmatrix

#include "k2tree.h"

typedef struct s_matrix {
   uint64_t elems; // how many filled cells
   uint logside; // log(max(width,height)
   uint transposed; // boolean
   uint64_t width,height; // matrix dimensions, filled w/0s up to 2^logside
   k2tree tree;
   } *matrix;

	// creates matrix of width x height with n cells (2n ints row,col) 
	// reorders cells array
matrix matCreate (uint64_t height, uint64_t width, uint64_t n, uint64_t *cells);

	// 32-bit version
matrix matCreate32 (uint64_t height, uint64_t width, uint64_t n, uint32_t *cells);

	// destroys matrix
void matDestroy (matrix M);
   
	// creates an empty matrix
matrix matEmpty (uint64_t height, uint64_t width);

	// creates a matrix with cell row,col
matrix matOne (uint64_t height, uint64_t width, uint64_t row, uint64_t col);

	// creates an identity matrix
matrix matId (uint64_t side);

        // creates a new copy of A, with its own data
matrix matCopy (matrix A);

	// transpose a matrix, creating a non-allocated copy that shares the 
	// data. You need not (and should not) matDestroy this copy
	// use *M = matTranspose(M) to actually transpose M
struct s_matrix matTranspose (matrix M);

        // writes M to file, which must be opened for writing
void matSave (matrix M, FILE *file);

        // loads matrix from file, which must be opened for reading
matrix matLoad (FILE *file);

	// space of the matrix, in w-bit words
uint64_t matSpace (matrix M);

	// dimensions of M, returns elems
uint64_t matDims (matrix M, uint *logside, uint64_t *width, uint64_t *height);

	// accesses a cell
uint matAccess (matrix M, uint64_t row, uint64_t col);

	// recovers all the cells in [r1..r2] x [c1..c2]
	// can use fullSide for r2 and c2
	// writes 2n integers in buffer, which must be preallocated
	// to 2*elems of uint64_t (worst case)
	// returns number of elements. just counts if buffer is NULL
uint64_t matCollect (matrix M, uint64_t r1, uint64_t r2, 
		     uint64_t c1, uint64_t c2, uint64_t *buffer);

	// 32-bit version
uint64_t matCollect32 (matrix M, uint32_t r1, uint32_t r2, 
		       uint32_t c1, uint32_t c2, uint32_t *buffer);

#define fullSide (~(uint64_t)0)
#define fullSide32 (~(uint)0)

        // (boolean) sum of two matrices, assumed to be of the same side
matrix matSum (matrix A, matrix B);

        // version with one row or one column, or both
matrix matSum1 (uint64_t row, matrix A, matrix B, uint64_t col);

// (boolean) sum of two matrices, assumed to be of the same side (recursive)
matrix matOr (matrix A, matrix B);

// version with one row or one column, or both (recursive)
matrix matOr1 (uint64_t row, matrix A, matrix B, uint64_t col);

        // (boolean) difference A-B, assumed to be of the same logside
matrix matDif (matrix A, matrix B);

        // version with one row or one column, or both
matrix matDif1 (uint64_t row, matrix A, matrix B, uint64_t col);

        // (boolean) symmetric difference, assumed to be of the same logside
matrix matXor (matrix A, matrix B);

        // version with one row or one column, or both
matrix matXor1 (uint64_t row, matrix A, matrix B, uint64_t col);

        // (boolean) intersection of A and B, assumed to be of the same logside
matrix matAnd (matrix A, matrix B);

        // version with one row or one column, or both
matrix matAnd1 (uint64_t row, matrix A, matrix B, uint64_t col);

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
matrix matMult (matrix A, matrix B);

	// version with one row or one column, or both
matrix matMult1 (uint64_t row, matrix A, matrix B, uint64_t col);

	// transitive closure of a matrix, pos says if it's + rather than *
matrix matClos (matrix A, uint pos);

	// versions to choose one row or one column, or both
matrix matClos1 (uint64_t row, matrix A, uint pos, uint64_t col);

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
matrix matMultClos1 (uint64_t row, matrix A, matrix B, uint pos, uint64_t col);

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
matrix matClosMult1 (uint64_t row, matrix A, uint pos, matrix B, uint64_t col);

        // multiplies M by vector V. allocates and returns the result M x V
        // use matrix transposition to do V^T x M = M^T x V (transposed vector)
        // V is assumed to be of size A->width, the output is of size A->height
double *matVectorMult (matrix A, double *vector);

#endif
