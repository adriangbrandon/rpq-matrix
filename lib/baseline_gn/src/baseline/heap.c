
#include "baseline/heap.h"

static void siftdown (heap H, uint n, uint i)

{ struct s_heap elem = H[i];
    while (1)
    { if (2*i+1 < n) // both children
        { if (H[2*i].key < H[2*i+1].key)
            { if (H[2*i].key < elem.key)
                { H[i] = H[2*i]; i = 2*i; }
                else break;
            }
            else
            { if (H[2*i+1].key < elem.key)
                { H[i] = H[2*i+1]; i = 2*i+1; }
                else break;
            }
        }
        else // left child or none
        { if ((2*i < n) && (H[2*i].key < elem.key))
            { H[i] = H[2*i]; i = 2*i; }
            break;
        }
    }
    H[i] = elem;
}

void heapify (heap H, uint n)

{ int i;
    for (i=n/2;i>=0;i--) siftdown (H,n,i);
}

extern inline struct s_heap findMin (heap H, uint n)

{ return H[0];
}

extern inline void replaceMin (heap H, uint n, uint key, uint data)

{ H[0].key = key;
    H[0].data = data;
    siftdown (H,n,0);
}

extern inline void extractMin (heap H, uint n)

{ n--;
    H[0] = H[n];
    siftdown (H,n,0);
}

