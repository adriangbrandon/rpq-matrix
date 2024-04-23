//
// Created by Adrián on 16/4/24.
//

#ifndef RPQ_MATRIX_BASELINE_HPP
#define RPQ_MATRIX_BASELINE_HPP

namespace bm_baseline {

    extern "C" {
    #include "baseline/matrix.h"
    #include "baseline/utilstime.h"
    }

    class wrapper {

    public:
        static const uint64_t full_side = fullSide;
        typedef matrix matrix_type;
        typedef s_matrix s_matrix_type;
        // creates matrix of width x height with n cells (2n ints row,col)
        // reorders cells array
        static inline matrix create(uint64_t height, uint64_t width, uint64_t n, uint *cells){
            return matCreate(height, width, n, cells);
        };

        static void time_begin(){
            time_beg();
        }

        // destroys matrix
        static inline void destroy(matrix M){
            matDestroy(M);
        };

        // creates an empty matrix
        static inline matrix empty(uint64_t height, uint64_t width){
            return matEmpty(height, width);
        };

        // creates an identity matrix
        static inline matrix id(uint64_t side){
            return matId(side);
        };

        // creates a new copy of A, with its own data
        static inline matrix copy(matrix A){
            return matCopy(A);
        };

        // transpose a matrix, creating a non-allocated copy that shares the
        // data. You need not (and should not) matDestroy this copy
        // use *M = matTranspose(M) to actually transpose M
        static inline struct s_matrix transpose(matrix M){
            return matTranspose(M);
        };

        // writes M to file, which must be opened for writing
        static inline void save(matrix M, FILE *file){
            matSave(M, file);
        };

        // loads matrix from file, which must be opened for reading
        static inline matrix load(FILE *file){
            return matLoad(file);
        };

        // space of the matrix, in w-bit words
        static inline uint64_t space(matrix M){
            return matSpace(M);
        };

        // dimensions of M, returns elems
        static inline uint64_t dims(matrix M, uint *logside, uint *width, uint *height){
            return matDims(M, width, height);
        };

        // accesses a cell
        static inline uint access(matrix M, uint64_t row, uint64_t col){
            return matAccess(M, row, col);
        }


        // (boolean) sum of two matrices, assumed to be of the same side
        static inline matrix sum (matrix A, matrix B){
            return matSum(A, B);
        };


        // version with one row or one column, or both
        static inline matrix sum1 (uint64_t row, matrix A, matrix B, uint64_t col){
            return matSum1(row, A, B, col);
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline matrix mult (matrix A, matrix B){
            return matMult(A, B);
        };

        // version with one row or one column, or both
        static inline matrix mult1 (uint64_t row, matrix A, matrix B, uint64_t col){
            return matMult1(row, A, B, col);
        };

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline matrix clos (matrix A, uint pos){
            return matClos(A, pos);
        };

        // versions to choose one row or one column, or both
        static inline matrix clos1 (uint64_t row, matrix A, uint pos, uint64_t col){
            return matClos1(row, A, pos, col);
        };

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
        static inline matrix mult_clos1 (uint64_t row, matrix A, matrix B, uint pos, uint64_t col){
            return matMultClos1(row, A, B, pos, col);
        };

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        static inline matrix clos_mult1 (uint64_t row, matrix A, uint pos, matrix B, uint64_t col){
            return matClosMult1(row, A, pos, B, col);
        };


    };
}


#endif //RPQ_MATRIX_BASELINE_HPP
