
/*

Nrgrep -- a fast and flexible pattern matching tool.
Copyright (C) 2000,2001 Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of 
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

*/

#ifndef REGULARINCLUDED
#define REGULARINCLUDED

/* Search a regular expression in a text */

#include "basics.hpp"  /* basics */
#include "options.hpp" /* buffer management */
#include "parser.hpp"  /* parser */
#include <string.h>

typedef struct

{
   int m; /* number of states, and of filtering one */

   int width;           /* shape of verif autom */
   int slices;           /* slices of  the verif autom */
   Mask initial, final; /* of fwd verif autom */

   Mask **fwdTrans; /* deterministic version of fwd verif automaton */
   Mask **bwdTrans; /* deterministic version of bwd verif automaton */
   Mask B[SIGMA];     /* B masks for verification */
   Mask V1;         /* extra space for verification */
} regularData;

/* Preprocesses pat and creates a regularData structure
	   for searching */

regularData *regularPreproc(const char *pat, Tree *tree, Tree **pos);
regularData *regularPreprocBInv(const char *pat, Tree *tree, Tree **pos, std::vector<int> &pos_pred);


/* Frees P */

void regularFree(regularData *P);

/* compute number of bits necessary to hold the regexp */

void regularLength(Tree *e, int *L);

/* reverses all the arrows of the NFA, exchanges initial and
           final states. the eps closures must have been done already */

void regularReverseArrows(Mask *trans, int m, Mask initial,
                          Mask final, Mask **rtrans, Mask *rinitial, Mask *rfinal);

/* builds a deterministic table for trans. it can handle multiwords
           in the det table too, but this is compatible with monoword usage.
           it receives the width of the deterministic table. the NFA is
           built so that the slices are always inside a single word */

void regularMakeDet(int width, Mask *trans, int m, Mask ***dtrans);

/* same as regularMakeDet but the masks are simple */

void regularMakeDet1(int width, Mask *trans, int m, mask ***dtrans);

/* receives pat and its strlen m, the number of bits L to store the
           regexp, trees e and pos[]. It fills the preallocated tables B,S,A,
           and the automaton: transitions *trans (not preallocated), and
           initial and final states (preallocated) */

void regularLoadMasks(const char *pat, int m, int L, Tree *e, Tree **pos,
                      Mask *B, Mask **trans, Mask initial, Mask final);

void regularLoadMasks(const char *pat, int m, int L, Tree *e, Tree **pos,
                      Mask *B, std::vector<int> &pos_pred,
                      Mask **trans, Mask initial, Mask final);

#endif
