
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

/* Search a regular expression in a text */

#include "../include/regular.hpp"

/* compute number of bits necessary to hold the regexp */

void regularLength(Tree *e, int *L)

{
	if (*L == 0)
		*L = 1;
	switch (e->type)
	{
	case STR:
		if (!e->eps)
			e->pos = (*L)++;
		break;
	case CONC:
	case OOR:
		regularLength(e->e1, L);
		regularLength(e->e2, L);
		break;
	case PLUS:
	case STAR:
	case QUESTION:
		regularLength(e->e1, L);
		break;
	}
}

static void firstLast(Tree *e, int L)

/* compute firstPos and lastPos, initial and final positions */

{
	switch (e->type)
	{
	case STR:
		e->firstPos = ZERO(createMask(L), L);
		if (!e->eps)
			e->firstPos = SET(e->firstPos, e->pos);
		e->lastPos = COPY(createMask(L), e->firstPos, L);
		break;
	case STAR:
	case PLUS:
	case QUESTION:
		firstLast(e->e1, L);
		e->firstPos = COPY(createMask(L), e->e1->firstPos, L);
		e->lastPos = COPY(createMask(L), e->e1->lastPos, L);
		break;
	case OOR:
		firstLast(e->e1, L);
		firstLast(e->e2, L);
		e->firstPos = OR(COPY(createMask(L), e->e1->firstPos, L),
						 e->e2->firstPos, L);
		e->lastPos = OR(COPY(createMask(L), e->e1->lastPos, L),
						e->e2->lastPos, L);
		break;
	case CONC:
		firstLast(e->e1, L);
		firstLast(e->e2, L);
		if (e->e1->eps)
			e->firstPos = OR(COPY(createMask(L), e->e1->firstPos, L),
							 e->e2->firstPos, L);
		else
			e->firstPos = COPY(createMask(L), e->e1->firstPos, L);
		if (e->e2->eps)
			e->lastPos = OR(COPY(createMask(L), e->e1->lastPos, L),
							e->e2->lastPos, L);
		else
			e->lastPos = COPY(createMask(L), e->e2->lastPos, L);
		break;
	}
}

static Mask follow(Tree *e, int pos, int L)

/* computes the follow set for a given position pos, that is,
	   the set of states which can follow position pos */

{
	Mask tmp1, tmp2;
	switch (e->type)
	{
	case STR:
		return ZERO(createMask(L), L);
		break;
	case STAR:
	case PLUS:
		if (ISSET(e->e1->lastPos, pos))
			return OR(follow(e->e1, pos, L), e->e1->firstPos, L);
		else
			return follow(e->e1, pos, L);
		break;
	case QUESTION:
		return follow(e->e1, pos, L);
		break;
	case OOR:
		if (ISSET(e->e1->maskPos, pos)) /* pos appears in left subtree */
		{
			if (ISSET(e->e2->maskPos, pos)) /* also in right subtree */
			{
				tmp1 = follow(e->e1, pos, L);
				tmp2 = follow(e->e2, pos, L);
				OR(tmp1, tmp2, L);
				free(tmp2);
				return tmp1;
			}
			else /* pos does not appear in right subtree */
				return follow(e->e1, pos, L);
		}
		else /* pos does not appear in left subtree */
		{
			if (ISSET(e->e2->maskPos, pos)) /* appears in right subtree */
				return follow(e->e2, pos, L);
			else /* does not appear anywhere */
				return ZERO(createMask(L), L);
		}
		break;
	case CONC:
		tmp1 = createMask(L);
		if (ISSET(e->e1->lastPos, pos))
			COPY(tmp1, e->e2->firstPos, L);
		else
			ZERO(tmp1, L);
		if (ISSET(e->e1->maskPos, pos)) /* pos appears in left subtree */
		{
			if (ISSET(e->e2->maskPos, pos)) /* also in right subtree */
			{
				tmp2 = follow(e->e1, pos, L);
				OR(tmp1, tmp2, L);
				free(tmp2);
				tmp2 = follow(e->e2, pos, L);
				OR(tmp1, tmp2, L);
				free(tmp2);
				return tmp1;
			}
			else /* pos does not appear in right subtree */
			{
				tmp2 = follow(e->e1, pos, L);
				OR(tmp1, tmp2, L);
				free(tmp2);
				return tmp1;
			}
		}
		else /* pos does not appear in left subtree */
		{
			if (ISSET(e->e2->maskPos, pos)) /* appears in right subtree */
			{
				tmp2 = follow(e->e2, pos, L);
				OR(tmp1, tmp2, L);
				free(tmp2);
				return tmp1;
			}
			else /* does not appear anywhere */
				return tmp1;
		}
		break;
	}
}

/* receives pat and its strlen m, the number of bits L to store the
	   regexp, trees e and pos[]. It fills the preallocated tables B,S,A,
	   and the automaton: transitions *trans (not preallocated), and
	   initial and final states (preallocated) */

void regularLoadMasks(const char *pat, int m, int L, Tree *e, Tree **pos,
					  Mask *B, Mask **trans, Mask initial, Mask final)

{
	int i, j;

	/* compute the B,S,A tables */
	i = 0;
	while (i < m)
		if (pos[i])
			getAclass(pat, &i, B, pos[i]->pos);
		else
			i++;
	if (OptRecChar != -1)
		ZERO(B[OptRecChar], L);

	/* compute maskPos */
	setMaskPos(e, L);
	/* compute firstPos and lastPos */
	firstLast(e, L);

	/* initial and final states */
	SET(ZERO(initial, L), 0);
	COPY(final, e->lastPos, L);

	/* If the expression admits the empty string then the initial states are also final */
	if (e->eps)
		OR(final, initial, L);

	/* allocate the transitions */
	*trans = (mask**)malloc(L * sizeof(Mask *));
	(*trans)[0] = createMasks(L, m);
	for (i = 1; i < L; i++)
	{
		(*trans)[i] = (*trans)[0] + i * maskSize(m);
		ZERO((*trans)[i], m);
	}

	/* load the transitions. this is basically the follow set, but it is
	    "first" for the initial state */
	for (i = 0; i < L; i++)
	{
		if (i == 0)
			COPY((*trans)[i], e->firstPos, L); /* initial position */
		else								   /* all other positions */
		{
			Mask tmp = follow(e, i, L);
			COPY((*trans)[i], tmp, L);
			free(tmp);
		}
	}
}

/* receives pat and its strlen m, the number of bits L to store the
	   regexp, trees e and pos[]. It fills the preallocated tables B,S,A,
	   and the automaton: transitions *trans (not preallocated), and
	   initial and final states (preallocated) */

void regularLoadMasks(const char *pat, int m, int L, Tree *e, Tree **pos,
                      Mask *B, std::vector<int> &pos_pred,
                      Mask **trans, Mask initial, Mask final)

{
    int i, j;

    /* compute the B,S,A tables */
    i = 0;
    while (i < m)
        if (pos[i])
            getAclass(pat, &i, B, pos[i]->pos, pos_pred);
        else
            i++;
    if (OptRecChar != -1)
        ZERO(B[OptRecChar], L);

    /* compute maskPos */
    setMaskPos(e, L);
    /* compute firstPos and lastPos */
    firstLast(e, L);

    /* initial and final states */
    SET(ZERO(initial, L), 0);
    COPY(final, e->lastPos, L);

    /* If the expression admits the empty string then the initial states are also final */
    if (e->eps)
        OR(final, initial, L);

    /* allocate the transitions */
    *trans = (mask**)malloc(L * sizeof(Mask *));
    (*trans)[0] = createMasks(L, m);
    for (i = 1; i < L; i++)
    {
        (*trans)[i] = (*trans)[0] + i * maskSize(m);
        ZERO((*trans)[i], m);
    }

    /* load the transitions. this is basically the follow set, but it is
        "first" for the initial state */
    for (i = 0; i < L; i++)
    {
        if (i == 0)
            COPY((*trans)[i], e->firstPos, L); /* initial position */
        else								   /* all other positions */
        {
            Mask tmp = follow(e, i, L);
            COPY((*trans)[i], tmp, L);
            free(tmp);
        }
    }
}

void regularReverseArrows(Mask *trans, int m, Mask initial,
						  Mask final, Mask **rtrans, Mask *rinitial, Mask *rfinal)

/* reverses all the arrows of the NFA, exchanges initial and
  	   final states. the eps closures must have been done already */

{
	int i, j;
	*rinitial = COPY(createMask(m), final, m);
	*rfinal = COPY(createMask(m), initial, m);
	*rtrans = (mask**)malloc(m * sizeof(Mask *));
	(*rtrans)[0] = createMasks(m, m);
	for (i = 0; i < m; i++)
	{
		(*rtrans)[i] = (*rtrans)[0] + i * maskSize(m);
		ZERO((*rtrans)[i], m);
	}
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			if (ISSET(trans[i], j))
				SET((*rtrans)[j], i);
}

void regularMakeDet(int width, Mask *trans, int m, Mask ***dtrans)

/* builds a deterministic table for trans. it can handle multiwords
	   in the det table too, but this is compatible with monoword usage.
	   it receives the width of the deterministic table. the NFA is
	   built so that the slices are always inside a single word */

{
	int slices;
	int i, j, f, w, b;
	mask m1, m2;
	int dm = 1 << width;
	slices = (m + width - 1) / width;
	/* allocate and structure the memory */
	*dtrans = (mask***)malloc(slices * sizeof(Mask **));
	(*dtrans)[0] = (mask**)malloc(slices * dm * sizeof(Mask *));
	for (i = 1; i < slices; i++)
		(*dtrans)[i] = (*dtrans)[0] + i * dm;
	(*dtrans)[0][0] = createMasks(slices * dm, m);
	for (i = 0; i < slices; i++)
	{
		(*dtrans)[i][0] = (*dtrans)[0][0] + i * dm * maskSize(m);
		for (j = 1; j < dm; j++)
			(*dtrans)[i][j] = (*dtrans)[i][0] + j * maskSize(m);
	}
	/* fill with NFA */
	f = 0;
	for (i = 0; i < slices; i++)
	{
		w = (i < slices - 1) ? width : m - f;
		ZERO((*dtrans)[i][0], m);
		for (b = 0; b < w; b++)
		{
			if (b + 1 < W)
			{
				m1 = ~(1 << b);
				m2 = 1 << (b + 1);
				for (j = 1 << b; j < m2; j++)
				{
					COPY((*dtrans)[i][j], (*dtrans)[i][j & m1], m);
					OR((*dtrans)[i][j], trans[f + b], m);
				}
			}
			else
			{
				m1 = ~(1 << b);
				for (j = 1 << b; j != 0; j++)
				{
					COPY((*dtrans)[i][j], (*dtrans)[i][j & m1], m);
					OR((*dtrans)[i][j], trans[f + b], m);
				}
			}
		}
		f += width;
	}
}

void regularMakeDet1(int width, Mask *trans, int m, mask ***dtrans)

/* same as regularMakeDet but the masks are simple */

{
	int slices;
	int i, j, f, w, b;
	mask m1, m2;
	int dm = 1 << width;
	slices = (m + width - 1) / width;
	/* allocate and structure the memory */
	*dtrans = (mask**)malloc(slices * sizeof(mask **));
	(*dtrans)[0] = (mask*)malloc(slices * dm * sizeof(mask));
	for (i = 1; i < slices; i++)
		(*dtrans)[i] = (*dtrans)[0] + i * dm;
	/* fill with NFA */
	f = 0;
	for (i = 0; i < slices; i++)
	{
		w = (i < slices - 1) ? width : m - f;
		(*dtrans)[i][0] = ZEROS;
		for (b = 0; b < w; b++)
		{
			if (b + 1 < W)
			{
				m1 = ~(1 << b);
				m2 = 1 << (b + 1);
				for (j = 1 << b; j < m2; j++)
					(*dtrans)[i][j] = (*dtrans)[i][j & m1] | *(trans[f + b]);
			}
			else
			{
				m1 = ~(1 << b);
				for (j = 1 << b; j != 0; j++)
					(*dtrans)[i][j] = (*dtrans)[i][j & m1] | *(trans[f + b]);
			}
		}
		f += width;
	}
}

regularData *regularPreproc(const char *pat, Tree *tree, Tree **pos)

{
    regularData *P = (regularData*)malloc(sizeof(regularData));
    int slices, i, j, c, *map;
    Mask *trans, *rtrans; /* nondet trans */

    /* allocate and load the masks */

    P->m = 0;
    regularLength(tree, &P->m);
    P->initial = createMask(P->m);
    P->final = createMask(P->m);
    P->B[0] = createMasks(SIGMA, P->m);
    for (c = 0; c < SIGMA; c++)
    {
        P->B[c] = ZERO(P->B[0] + c * maskSize(P->m), P->m);
    }

    regularLoadMasks(pat, strlen(pat), P->m, tree, pos, P->B, &trans, P->initial, P->final);

    /* make rtrans = reverse of trans, O(m^2/w + m^2) time */

    rtrans = (mask**) malloc(P->m * sizeof(Mask));
    rtrans[0] = createMasks(P->m, P->m);
    for (i = 0; i < P->m; i++)
        rtrans[i] = ZERO(rtrans[0] + i * maskSize(P->m), P->m);
    for (i = 0; i < P->m; i++)
        for (j = 0; j < P->m; j++)
            if (ISSET(trans[i], j))
                SET(rtrans[j], i);

    /* create verification automata */
    slices = (P->m + OptDetWidth - 1) / OptDetWidth;
    P->slices = slices;
    P->width = (P->m + slices - 1) / slices;
    regularMakeDet(P->width, trans, P->m, &P->fwdTrans);
    regularMakeDet(P->width, rtrans, P->m, &P->bwdTrans);

    free(trans[0]);
    free(trans);
    free(rtrans[0]);
    free(rtrans);
    P->V1 = createMask(P->m);
    return P;
}


regularData *regularPreprocBInv(const char *pat, Tree *tree, Tree **pos, std::vector<int> &pos_pred)

{
	regularData *P = (regularData*)malloc(sizeof(regularData));
	int slices, i, j, c, *map;
	Mask *trans, *rtrans; /* nondet trans */

	/* allocate and load the masks */

	P->m = 0;
	regularLength(tree, &P->m);
	P->initial = createMask(P->m);
	P->final = createMask(P->m);
	P->B[0] = createMasks(SIGMA, P->m);
	for (c = 0; c < SIGMA; c++)
	{
		P->B[c] = ZERO(P->B[0] + c * maskSize(P->m), P->m);
	}

	regularLoadMasks(pat, strlen(pat), P->m, tree, pos, P->B, pos_pred, &trans, P->initial, P->final);

	/* make rtrans = reverse of trans, O(m^2/w + m^2) time */

	rtrans = (mask**) malloc(P->m * sizeof(Mask));
	rtrans[0] = createMasks(P->m, P->m);
	for (i = 0; i < P->m; i++)
		rtrans[i] = ZERO(rtrans[0] + i * maskSize(P->m), P->m);
	for (i = 0; i < P->m; i++)
		for (j = 0; j < P->m; j++)
			if (ISSET(trans[i], j))
				SET(rtrans[j], i);

	/* create verification automata */
	slices = (P->m + OptDetWidth - 1) / OptDetWidth;
	P->slices = slices;
	P->width = (P->m + slices - 1) / slices;
	regularMakeDet(P->width, trans, P->m, &P->fwdTrans);
	regularMakeDet(P->width, rtrans, P->m, &P->bwdTrans);

	free(trans[0]);
	free(trans);
	free(rtrans[0]);
	free(rtrans);
	P->V1 = createMask(P->m);
 	return P;
}

void regularFree(regularData *P)

/* Frees P */

{
	free(P->B[0]);
	free(P->fwdTrans[0][0]);
	free(P->fwdTrans[0]);
	free(P->fwdTrans);
	free(P->initial);
	free(P->final);
	free(P->bwdTrans[0][0]);
	free(P->bwdTrans[0]);
	free(P->bwdTrans);
	free(P->V1);
	free(P);
}
