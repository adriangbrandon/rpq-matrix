
	// simple linear probing hash table, using twice the needed space
	// ad-hoc: the values are the positions of creation array

#include "basics.h"

#define noval (~(uint)0)

typedef struct {
   uint *table;
   uint *sat;
   uint b; // mask to map to table, (1 << smth)-1
   } hashTable;

	// searches for key in hash and returns value, or noval if not found

extern inline uint hashSearch (hashTable hash, uint key);

	// inserts (key,val) in hash table, if it's not there

extern inline void hashInsert (hashTable hash, uint key, uint val);

	// creates and returns hash table for distinct values[0..n-1]
	// the value of value[i] is i

hashTable hashCreate (uint *values, uint n);

	// destroys hash

void hashDestroy (hashTable hash);

