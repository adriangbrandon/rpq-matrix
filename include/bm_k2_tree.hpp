//
// Created by Adrián on 16/4/24.
//

#ifndef RPQ_MATRIX_K2_TREE_HPP
#define RPQ_MATRIX_K2_TREE_HPP

namespace bm_k2_tree {

    extern "C" {
    #include "k2-tree/matrix.h"
    #include "k2-tree/utilstime.h"
    }

    class wrapper {

    public:
        static const uint64_t full_side = fullSide;
        typedef matrix matrix_type;
        typedef s_matrix s_matrix_type;
        // creates matrix of width x height with n cells (2n ints row,col)
        // reorders cells array
        static inline matrix_type create(uint64_t height, uint64_t width, uint64_t n, uint64_t *cells){
            return matCreate(height, width, n, cells);
        };

        static void time_begin(){
            time_beg();
        }

        // 32-bit version
        static inline matrix_type create32 (uint64_t height, uint64_t width, uint64_t n, uint32_t *cells){
            return matCreate32(height, width, n, cells);
        }

        // destroys matrix
        static inline void destroy(matrix_type M){
            matDestroy(M);
        };

        // creates an empty matrix
        static inline matrix_type empty(uint64_t height, uint64_t width){
            return matEmpty(height, width);
        };

        // creates a matrix with cell row,col
        static inline matrix_type one(uint64_t height, uint64_t width, uint64_t row, uint64_t col){
            return matOne(height, width, row, col);
        };

        // creates an identity matrix
        static inline matrix_type id(uint64_t side){
            return matId(side);
        };

        // creates a new copy of A, with its own data
        static inline matrix_type copy(matrix A){
            return matCopy(A);
        };

        // transpose a matrix, creating a non-allocated copy that shares the
        // data. You need not (and should not) matDestroy this copy
        // use *M = matTranspose(M) to actually transpose M
        static inline s_matrix_type transpose(matrix_type M){
            return matTranspose(M);
        };

        // writes M to file, which must be opened for writing
        static inline void save(matrix_type M, FILE *file){
            matSave(M, file);
        };

        // loads matrix from file, which must be opened for reading
        static inline matrix_type load(FILE *file){
            return matLoad(file);
        };

        // space of the matrix, in w-bit words
        static inline uint64_t space(matrix_type M){
            return matSpace(M);
        };

        // dimensions of M, returns elems
        static inline uint64_t dims(matrix_type M, uint *logside, uint64_t *width, uint64_t *height){
            return matDims(M, logside, width, height);
        };

        // accesses a cell
        static inline uint access(matrix_type M, uint64_t row, uint64_t col){
            return matAccess(M, row, col);
        }

        // recovers all the cells in [r1..r2] x [c1..c2]
        // can use fullSide for r2 and c2
        // writes 2n integers in buffer, which must be preallocated
        // to 2*elems of uint64_t (worst case)
        // returns number of elements. just counts if buffer is NULL
        static inline uint64_t collect (matrix_type M, uint64_t r1, uint64_t r2,
                                    uint64_t c1, uint64_t c2, uint64_t *buffer){
            return matCollect(M, r1, r2, c1, c2, buffer);
        };

        // 32-bit version
        static inline uint64_t collect32 (matrix_type M, uint32_t r1, uint32_t r2,
                                      uint32_t c1, uint32_t c2, uint32_t *buffer){
            return matCollect32(M, r1, r2, c1, c2, buffer);
        };


        // (boolean) sum of two matrices, assumed to be of the same side
        static inline matrix_type sum (matrix_type A, matrix_type B){
            return matSum(A, B);
        };

        static inline matrix_type mat_or (matrix_type A, matrix_type B){
            return matOr(A, B);
        }; // recursive implementation

        // version with one row or one column, or both
        static inline matrix_type sum1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matSum1(row, A, B, col);
        };
        static inline matrix_type mat_or1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matOr1(row, A, B, col);
        }; // recursive

        // (boolean) difference A-B, assumed to be of the same logside
        static inline matrix_type dif (matrix_type A, matrix_type B){
            return matDif(A, B);
        };

        // version with one row or one column, or both
        static inline matrix_type dif1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matDif1(row, A, B, col);
        }

        // (boolean) symmetric difference, assumed to be of the same logside
        static inline matrix_type mat_xor (matrix A, matrix B){
            return matXor(A, B);
        };

        // version with one row or one column, or both
        static inline matrix_type mat_xor1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matXor1(row, A, B, col);
        };

        // (boolean) intersection of A and B, assumed to be of the same logside
        static inline matrix_type mat_and (matrix_type A, matrix_type B){
            return matAnd(A, B);
        };

        // version with one row or one column, or both
        static inline matrix_type mat_and1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matAnd1(row, A, B, col);
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline matrix_type mult (matrix_type A, matrix_type B){
            return matMult(A, B);
        };

        // version with one row or one column, or both
        static inline matrix_type mult1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matMult1(row, A, B, col);
        };

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline matrix_type clos (matrix_type A, uint pos){
            return matClos(A, pos);
        };

        // versions to choose one row or one column, or both
        static inline matrix_type clos1 (uint64_t row, matrix_type A, uint pos, uint64_t col){
            return matClos1(row, A, pos, col);
        };

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
        static inline matrix_type mult_clos1 (uint64_t row, matrix_type A, matrix_type B, uint pos, uint64_t col){
            return matMultClos1(row, A, B, pos, col);
        };

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        static inline matrix_type clos_mult1 (uint64_t row, matrix_type A, uint pos, matrix_type B, uint64_t col){
            return matClosMult1(row, A, pos, B, col);
        };


    };
}


#endif //RPQ_MATRIX_K2_TREE_HPP
