//
// Created by Adrián on 2/5/23.
//

// reads the pairs/*.pairs files and creates *.mat files with their
// matrix representations

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

extern "C" {
    #include <matrix.h>
}

#define N 958844164
#define S 5420 // 1 to 5419
#define V 296008192 // 1 to...

// #define MAIN


uint64_t cellsA[] = { 1,2, 2,2, 1,0, 0,3, 5,2, 6,2, 5,0, 4,3,
                      1,6, 2,6, 1,4, 0,7, 5,6, 6,6, 5,4, 4,7 };

uint64_t cellsB[] = { 3,3, 1,0, 0,3, 6,2, 5,0, 7,7 };

int main (void)

{ matrix A = matCreate(8,8,16,cellsA);
    matrix B = matCreate(8,8,6,cellsB);
    matrix Ss,M,C;
    uint64_t row,col;
    uint64_t elems,width,height;
    uint logside;

    elems = matDims (A,&logside,&width,&height);
    printf("Matrix A: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(A,row,col));
        printf("\n");
    }
    printf ("\n");
    elems = matDims (B,&logside,&width,&height);
    printf("Matrix B: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(B,row,col));
        printf("\n");
    }
    printf ("\n");
    matTranspose(B); printf("B transposed\n");
    elems = matDims (B,&logside,&width,&height);
    printf("Matrix B: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(B,row,col));
        printf("\n");
    }
    printf ("\n");
    Ss = matSum(A,B);
    elems = matDims (Ss,&logside,&width,&height);
    printf("Boolean sum: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(Ss,row,col));
        printf("\n");
    }
    printf ("\n");
//    Ss = matSum1(A,B,3,3);
//    elems = matDims (Ss,&logside,&width,&height);
//    printf("Boolean sum at 3,3: %lix%li, %i bits, %li elements\n\n",
//	   height,width,logside,elems);
//    for (row=0;row<8;row++)
//        { for (col=0;col<8;col++)
//	      printf("%i",matAccess(Ss,row,col));
//	  printf("\n");
//	}
//    printf ("\n");
    M = matMult(A,B);
    elems = matDims (M,&logside,&width,&height);
    printf("Boolean product: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { int k;
        for (col=0;col<8;col++)
            printf("%i",matAccess(M,row,col));
        printf("    ");
        for (col=0;col<8;col++)
        { int v = 0;
            for (k=0;k<8;k++)
                v |= matAccess(A,row,k) & matAccess(B,k,col);
            printf ("%i",v);
        }
        printf("\n");
    }
    printf ("\n");
    C = matClos(Ss,1);
    elems = matDims (C,&logside,&width,&height);
    printf("Closure of sum: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(C,row,col));
        printf("\n");
    }
    printf ("\n");
    C = matClos1(Ss,1,fullSide,4);
    elems = matDims (C,&logside,&width,&height);
    printf("Closure of sum, col 5: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(C,row,col));
        printf("\n");
    }
    printf ("\n");
    C = matClos1(Ss,1,1,fullSide);
    elems = matDims (C,&logside,&width,&height);
    printf("Closure of sum, row 2: %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(C,row,col));
        printf("\n");
    }
    printf ("\n");
    C = matClos1(Ss,1,1,4);
    elems = matDims (C,&logside,&width,&height);
    printf("Closure of sum, cell (2,5): %lix%li, %i bits, %li elements\n\n",
           height,width,logside,elems);
    for (row=0;row<8;row++)
    { for (col=0;col<8;col++)
            printf("%i",matAccess(C,row,col));
        printf("\n");
    }
    printf ("\n");
}
