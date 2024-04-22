
// min-heap data structure, rough

#include "basics.h"

typedef struct s_heap
{ uint key;
    uint data;
} *heap;

// takes array H and gives it heap order

void heapify (heap H, uint n);

// returns the key

extern inline struct s_heap findMin (heap H, uint n);

// extracts min, n must decrease

extern inline void extractMin (heap H, uint n);

// replaces and sifts down

extern inline void replaceMin (heap H, uint n, uint key, uint data);
