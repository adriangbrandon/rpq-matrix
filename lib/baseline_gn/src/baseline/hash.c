
#include "baseline/hash.h"

#define prime 4294967291 // 32-bit unsigned prime

extern inline uint hashSearch (hashTable hash, uint key)

   { uint pos = (key * prime) & hash.b;
     uint tkey;
     while ((tkey = hash.table[pos]) != key)
        { if (tkey == ~0) return noval; // not found;
          pos = (pos+1) & hash.b;
        }
     return hash.sat[pos];
   }

extern inline void hashInsert (hashTable hash, uint key, uint val)

   { uint pos = (key * prime) & hash.b;
     while (hash.table[pos] != ~0)
        pos = (pos+1) & hash.b;
     hash.table[pos] = key;
     hash.sat[pos] = val;
   }

hashTable hashCreate (uint *values, uint n)

   { hashTable hash;
    uint size = 1 << (1 + numbits(n));
    hash.b = size - 1;
    hash.table = (uint*)myalloc(size * sizeof(uint));
    hash.sat   = (uint*)myalloc(size * sizeof(uint));
    for (uint i = 0; i < size; i++) hash.table[i] = ~0;
    for (uint i = 0; i < n; i++) hashInsert(hash, values[i], i);
    return hash;
   }

void hashDestroy (hashTable hash)

   { myfree(hash.table);
     myfree(hash.sat);
   }

