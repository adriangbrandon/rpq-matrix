//
// Created by Adrián on 16/4/24.
//

#ifndef RPQ_MATRIX_BASELINE_32_HPP
#define RPQ_MATRIX_BASELINE_32_HPP

namespace bm_baseline_32 {

    extern "C" {
    #include "baseline/matrix32.h"
    #include "baseline/utilstime.h"
    }

    class wrapper {

    public:

        static const uint64_t full_side = fullSide;
        typedef matrix32 matrix_type;
        typedef s_matrix32 s_matrix_type;
        // creates matrix of width x height with n cells (2n ints row,col)
        // reorders cells array
        static inline matrix_type create(uint64_t height, uint64_t width, uint64_t n, uint *cells){
            return matCreate32(height, width, n, cells);
        };

        static void time_begin(){
            time_beg();
        }

        // destroys matrix
        static inline void destroy(matrix_type M){
            matDestroy32(M);
        };

        // creates an empty matrix
        static inline matrix_type empty(uint64_t height, uint64_t width){
            return matEmpty32(height, width);
        };

        // creates an identity matrix
        static inline matrix_type id(uint64_t side){
            return matId32(side);
        };

        // version of identity matrix with one row or column
        static inline matrix_type id1(uint64_t side, uint64_t rc){
            return matOne32(side, side, rc, rc);
        };

        // creates a new copy of A, with its own data
        static inline matrix_type copy(matrix_type A){
            return matCopy32(A);
        };

        // transpose a matrix, creating a non-allocated copy that shares the
        // data. You need not (and should not) matDestroy this copy
        // use *M = matTranspose(M) to actually transpose M
        static inline s_matrix_type transpose(matrix_type M){
            return matTranspose32(M);
        };

        // writes M to file, which must be opened for writing
        static inline void save(matrix_type M, FILE *file){
            matSave32(M, file);
        };

        // loads matrix from file, which must be opened for reading
        static inline matrix_type load(FILE *file){
            return matLoad32(file);
        };

        // space of the matrix, in w-bit words
        static inline uint64_t space(matrix_type M){
            return matSpace32(M);
        };

        // dimensions of M, returns elems
        static inline uint64_t dims(matrix_type M, uint *logside, uint *width, uint *height){
            return matDims32(M, width, height);
        };

        // accesses a cell
        static inline uint access(matrix_type M, uint64_t row, uint64_t col){
            return matAccess32(M, row, col);
        }


        // (boolean) sum of two matrices, assumed to be of the same side
        static inline matrix_type sum (matrix_type A, matrix_type B){
            return matSum32(A, B);
        };


        // version with one row or one column, or both
        static inline matrix_type sum1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matSum132(row, A, B, col);
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline matrix_type mult (matrix_type A, matrix_type B){
            return matMult32(A, B);
        };

        // version with one row or one column, or both
        static inline matrix_type mult1 (uint64_t row, matrix_type A, matrix_type B, uint64_t col){
            return matMult132(row, A, B, col);
        };

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline matrix_type clos (matrix_type A, uint pos){
            return matClos32(A, pos);
        };

        // versions to choose one row or one column, or both
        static inline matrix_type clos1 (uint64_t row, matrix_type A, uint pos, uint64_t col){
            return matClos132(row, A, pos, col);
        };

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
        static inline matrix_type mult_clos1 (uint64_t row, matrix_type A, matrix_type B, uint pos, uint64_t col){
            return matMultClos132(row, A, B, pos, col);
        };

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        static inline matrix_type clos_mult1 (uint64_t row, matrix_type A, uint pos, matrix_type B, uint64_t col){
            return matClosMult132(row, A, pos, B, col);
        };


    };
}


#endif //RPQ_MATRIX_BASELINE_32_HPP
