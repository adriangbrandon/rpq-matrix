
#ifndef INCLUDEDk2tree
#define INCLUDEDk2tree

#include "bitvector.h"

typedef struct s_k2tree {
   uint nlevels;
   bitvector B;
   } *k2tree;

typedef uint64_t k2node;

#define k2null (~0)

	// creates k2tree from n coords (2n ints x,y) that need nbits bits
	// reorders coords array
k2tree k2create (uint64_t n, uint nbits, uint64_t *coords);

	// 32-bit version
k2tree k2create32 (uint64_t n, uint nbits, uint32_t *coords);

	// creates k2tree from bits, of length len, nlevels levels
	// if own, will free bitvector when destroyed
k2tree k2createFrom (uint nlevels, uint64_t len, void *bits, uint own);

	// destroys k2tree, frees bitvector if owned
void k2destroy (k2tree T);

        // creates a new copy of T, with its own data
k2tree k2copy (k2tree T);

        // writes T to file, which must be opened for writing
void k2save (k2tree T, FILE *file);

        // loads k2tree from file, which must be opened for reading
k2tree k2load (FILE *file);

	// space of k2tree in w-bit words
uint64_t k2space (k2tree T);

	// levels of the k2tree
extern inline uint k2levels (k2tree T);

	// bitvector of the k2tree (somehow internal)
extern inline bitvector k2bits (k2tree T);

	// root of the k2tree
extern inline k2node k2root (k2tree T);

	// u has i-th child, assumes u is right and i = 0..3
extern inline int k2hasChild (k2tree T, k2node u, uint i);

	// i-th child of u, assumes u is right, i = 0..3, and child exists
extern inline int k2child (k2tree T, k2node u, uint i);

        // recovers all the cells in [r1..r2] x [c1..c2]
        // writes 2n integers in buffer, which must be preallocated
        // to 2*elems of uint64_t (worst case)
	// returns number of elements found. just counts if buffer is NULL
	// cr says to write col,row, otherwise it is row,col
uint64_t k2collect (k2tree T, uint64_t r1, uint64_t r2, 
		    uint64_t c1, uint64_t c2, uint64_t *buffer, uint cr);

	// 32-bit version
uint64_t k2collect32 (k2tree T, uint32_t r1, uint32_t r2, 
		    uint32_t c1, uint32_t c2, uint32_t *buffer, uint cr);

        // merges the bit arrays of two k2trees, writes its bit length
        // in *len and number of elements in *telems. nodes per level are
        // written in array *levels if not null, levels[0] is #leaves
	// Btransp = 1 if B has to be interpreted as transposed
uint64_t *k2merge (uint64_t *treeA, uint64_t lenA,
                   uint64_t *treeB, uint64_t lenB, uint level, 
		   uint64_t *len, uint64_t *telems, uint64_t *levels);

	// same but B is interpreted as transposed
uint64_t *k2mergeT (uint64_t *treeA, uint64_t lenA, k2tree B,
		    uint64_t *len, uint64_t *telems, uint64_t *levels);

#endif
