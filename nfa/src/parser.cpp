
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

#include "../include/parser.hpp"
#include "../include/options.hpp"
#include <iostream>

#define Sep '#'
#define All '.'
#define OpenPar '('
#define ClosePar ')'
#define Star '*'
#define Plus '+'
#define Question '?'
#define Or '|'
#define OpenBracket '['
#define CloseBracket ']'
#define BackSlash '\\'
#define Dash '-'
#define Neg '^'

int getAchar(const char *pat, int *i)

/* Gets one character assuming it is not special.
	   Returns -1 in case of syntax error: \x not followed by 2 hex
	   digits or \ terminating the string. */
/* Syntax: \n -> newline, \t -> tab, \xXY -> hex ascii code XY,
		   \X -> literal character X, else normal char */

{
	int c;
	if (pat[*i] == BackSlash)
	{
		(*i)++;
		if (pat[*i] == 'n')
			c = '\n';
		else if (pat[*i] == 't')
			c = '\t';
		else if (tolower(pat[*i]) == 'x')
		{
			(*i)++;
			if (isdigit(pat[*i]))
				c = pat[*i] - '0';
			else if (isxdigit(pat[*i]))
				c = (tolower(pat[*i]) - 'a') + 10;
			else
				return -1;
			(*i)++;
			c <<= 4;
			if (isdigit(pat[*i]))
				c += pat[*i] - '0';
			else if (isxdigit(pat[*i]))
				c += (tolower(pat[*i]) - 'a') + 10;
			else
				return -1;
		}
		else if (pat[*i])
			c = pat[*i];
		else
			return -1;
	}
	else
		c = pat[*i];
	(*i)++;
	return c;
}

boolean getAclass(const char *pat, int *i, Mask *B, int l)

/* Gets a class of characters from pat[i...] (sets i at the next
	   position to read) and sets the appropriate bits in
	   B[*] column l. It assumes that the first character is not a
	   metacharacter. It returns true if it could read a class, false
	   in case of syntax error: from getAchar or missing ]. It can
	   also be used to test syntax, with B = NULL */
/* Syntax: class -> [set*] or [^set*] (^ complements the set)
		   set -> X-Y  (range of chars) or X (single char)
		   X,Y -> . (a class with all), # (a separator) or getAchar
	*/

{
	int c, c1, c2;
	boolean pos = true;
	switch (pat[*i])
	{
	case OpenBracket: /* a class is opened */
		(*i)++;
		if (pat[*i] == Neg) /* reverse meaning */
		{
			if (B)
				for (c = 0; c < SIGMA; c++)
					SET(B[c], l);
			pos = false;
			(*i)++;
		}
		while (pat[*i] && (pat[*i] != CloseBracket))
		{
			c1 = getAchar(pat, i);
			if (c1 == -1)
				return false;
			if ((pat[*i] == Dash) && pat[*i + 1] &&
				(pat[*i + 1] != CloseBracket))
			{
				(*i)++;
				c2 = getAchar(pat, i);
				if (c2 == -1)
					return false;
			}
			else
				c2 = c1;
			if (B)
				for (c = c1; c <= c2; c++)
				{
					if (pos)
						SET(B[c], l);
					else
						CLEAR(B[c], l);
				}
		}
		if (!pat[*i])
			return false; /* ] missing */
		(*i)++;
		break;
	case All: /* a class with everything */
		if (B)
			for (c = 0; c < SIGMA; c++)
				SET(B[c], l);
		(*i)++;
		break;
	case Sep: /* a class of separators */
		if (B)
			for (c = 0; c < SIGMA; c++)
				if (!isalnum(c))
					SET(B[c], l);
		(*i)++;
		break;
	case OpenPar:
	case ClosePar:
	case Or:
	case Question:
	case Star:
	case Plus:
		/* higher level constructions, reject */
		return false;
		break;
	default: /* a plain letter or escape code */
		c = getAchar(pat, i);
		if (c == -1)
			return false;
		if (B)
			SET(B[c], l);
		break;
	}

	/* expand for case insensitive searching */
	if (B && OptCaseInsensitive)
	{
		for (c = 'A'; c <= 'Z'; c++)
			if (pos && ISSET(B[c], l))
				SET(B[tolower(c)], l);
			else if (!pos && !ISSET(B[c], l))
				CLEAR(B[tolower(c)], l);
		for (c = 'a'; c <= 'z'; c++)
			if (pos && ISSET(B[c], l))
				SET(B[toupper(c)], l);
			else if (!pos && !ISSET(B[c], l))
				CLEAR(B[toupper(c)], l);
	}

	return true;
}

boolean getAclass(const char *pat, int *i, Mask *B, int l, int* pos_pred)

/* Gets a class of characters from pat[i...] (sets i at the next
	   position to read) and sets the appropriate bits in
	   B[*] column l. It assumes that the first character is not a
	   metacharacter. It returns true if it could read a class, false
	   in case of syntax error: from getAchar or missing ]. It can
	   also be used to test syntax, with B = NULL */
/* Syntax: class -> [set*] or [^set*] (^ complements the set)
		   set -> X-Y  (range of chars) or X (single char)
		   X,Y -> . (a class with all), # (a separator) or getAchar
	*/

{
	int c, c1, c2;
	boolean pos = true;
	switch (pat[*i])
	{
		case OpenBracket: /* a class is opened */
			(*i)++;
			if (pat[*i] == Neg) /* reverse meaning */
			{
				if (B)
					for (c = 0; c < SIGMA; c++) {
						SET(B[c], l);
						pos_pred[l] = c;
					}
				pos = false;
				(*i)++;
			}
			while (pat[*i] && (pat[*i] != CloseBracket))
			{
				c1 = getAchar(pat, i);
				if (c1 == -1)
					return false;
				if ((pat[*i] == Dash) && pat[*i + 1] &&
					(pat[*i + 1] != CloseBracket))
				{
					(*i)++;
					c2 = getAchar(pat, i);
					if (c2 == -1)
						return false;
				}
				else
					c2 = c1;
				if (B)
					for (c = c1; c <= c2; c++)
					{
						if (pos){
							SET(B[c], l);
							pos_pred[l] = c;
						}
						else {
							CLEAR(B[c], l);
							pos_pred[l] = 0;
						}
					}
			}
			if (!pat[*i])
				return false; /* ] missing */
			(*i)++;
			break;
		case All: /* a class with everything */
			if (B)
				for (c = 0; c < SIGMA; c++){
					SET(B[c], l);
					pos_pred[l] = c;
				}
			(*i)++;
			break;
		case Sep: /* a class of separators */
			if (B)
				for (c = 0; c < SIGMA; c++) {
					if (!isalnum(c)) {
						SET(B[c], l);
						pos_pred[l] = c;
					}
				}
			(*i)++;
			break;
		case OpenPar:
		case ClosePar:
		case Or:
		case Question:
		case Star:
		case Plus:
			/* higher level constructions, reject */
			return false;
			break;
		default: /* a plain letter or escape code */
			c = getAchar(pat, i);
			if (c == -1)
				return false;
			if (B) {
				SET(B[c], l);
				pos_pred[l] = c;
			}
			break;
	}

	/* expand for case insensitive searching */
	if (B && OptCaseInsensitive)
	{
		for (c = 'A'; c <= 'Z'; c++)
			if (pos && ISSET(B[c], l)) {
				SET(B[tolower(c)], l);
				pos_pred[l] = tolower(c);
			}
			else if (!pos && !ISSET(B[c], l)) {
				CLEAR(B[tolower(c)], l);
				pos_pred[l] = 0;
			}
		for (c = 'a'; c <= 'z'; c++)
			if (pos && ISSET(B[c], l)) {
				SET(B[toupper(c)], l);
				pos_pred[l] = toupper(c);
			}
			else if (!pos && !ISSET(B[c], l)) {
				CLEAR(B[toupper(c)], l);
				pos_pred[l] = 0;
			}
	}

	return true;
}

/* Regexp Syntax:
		inside a class [...]
			^ at the beginning: negation of class
			\c: take character c as literal
			\n,\t: newline, tab
			\xdd: hex ascii code
	        outside a class, in addition
	                <class>+
	                <class>?
			<class>*
			# = separator
	                . = a class of all                                      
		unique of regexps
			(exp), exp|exp, exp+, exp*, exp?, exp exp
	*/

void freeTree(Tree *e)

/* frees a structure allocated by parse */

{
	if (e == NULL)
		return;
	free(e->maskPos);
	free(e->firstPos);
	free(e->lastPos);
	freeTree(e->e1);
	freeTree(e->e2);
	free(e);
}

static Tree *parseOr(const char *pat, int *i, Tree **pos);

static Tree *parseOne(const char *pat, int *i, Tree **pos)

/* parses one subexpression */

{
	Tree *e;
	if (pat[*i] == OpenPar) /* opens subexpression */
	{
		(*i)++;
		e = parseOr(pat, i, pos);
		if (e == NULL)
			return NULL;
		if (pat[*i] != ClosePar)
		{
			freeTree(e);
			return NULL;
		}
		(*i)++;
	}
	else
	{
		e = (Tree*)malloc(sizeof(Tree));
		e->type = STR;
		e->e1 = e->e2 = NULL;
		/* could be the empty string if higher level chars are here */
		e->eps = !pat[*i] || (pat[*i] == ClosePar) || (pat[*i] == Star) ||
				 (pat[*i] == Plus) || (pat[*i] == Question) || (pat[*i] == Or);
		if (!e->eps)
		{
			pos[*i] = e;
			if (!getAclass(pat, i, NULL, 0))
			{
				free(e);
				return NULL;
			}
		}
		e->maskPos = e->firstPos = e->lastPos = NULL;
	}
	return e;
}


static Tree *parseOneClosed(const char *pat, int *i, Tree **pos)

/* parses one subexpression plus closures +,*,? */

{
	Tree *e, *e1;

	e = parseOne(pat, i, pos);
	if (e == NULL)
		return NULL;

	while (true)
		switch (pat[*i])
		{
		case Plus: /* plus */
			e1 = (Tree*)malloc(sizeof(Tree));
			e1->type = PLUS;
			e1->eps = e->eps;
			e1->e1 = e;
			e1->e2 = NULL;
			e = e1;
			e->maskPos = e->firstPos = e->lastPos = NULL;
			(*i)++;
			break;
		case Star: /* star */
			e1 = (Tree*)malloc(sizeof(Tree));
			e1->type = STAR;
			e1->eps = true;
			e1->e1 = e;
			e1->e2 = NULL;
			e = e1;
			e->maskPos = e->firstPos = e->lastPos = NULL;
			(*i)++;
			break;
		case Question: /* previous is optional */
			e1 = (Tree*)malloc(sizeof(Tree));
			e1->type = QUESTION;
			e1->eps = true;
			e1->e1 = e;
			e1->e2 = NULL;
			e = e1;
			e->maskPos = e->firstPos = e->lastPos = NULL;
			(*i)++;
			break;
		default:
			return e; /* end of closed expression */
		}
}



static Tree *parseConc(const char *pat, int *i, Tree **pos)

/* parses a concatenation */

{
	Tree *e, *e1;
	e = parseOneClosed(pat, i, pos);
	if (e == NULL)
		return NULL;
	switch (pat[*i])
	{
	case Or:
	case ClosePar:
	case 0: /*or with next or end of pattern */
		break;
	default: /* must be a concatenation */
		e1 = (Tree*)malloc(sizeof(Tree));
		e1->type = CONC;
		e1->e1 = e;
		e = e1;
		e->e2 = parseConc(pat, i, pos);
		if (e->e2 == NULL)
		{
			freeTree(e);
			return NULL;
		}
		e->eps = e->e1->eps && e->e2->eps;
		e->maskPos = e->firstPos = e->lastPos = NULL;
		break;
	}
	return e;
}



static Tree *parseOr(const char *pat, int *i, Tree **pos)

/* parses a disjunction */

{
	Tree *e, *e1;
	e = parseConc(pat, i, pos);
	if (e == NULL)
		return NULL;
	switch (pat[*i])
	{
	case ClosePar:
	case 0: /* end of pattern */
		break;
	default: /* must be an or */
		e1 = (Tree*)malloc(sizeof(Tree));
		e1->type = OOR;
		e1->e1 = e;
		e = e1;
		(*i)++; /* skip the | */
		e->e2 = parseOr(pat, i, pos);
		if (e->e2 == NULL)
		{
			freeTree(e);
			return NULL;
		}
		e->eps = e->e1->eps || e->e2->eps;
		e->maskPos = e->firstPos = e->lastPos = NULL;
		break;
	}
	return e;
}


static void simpFree(Tree *e, Tree **pos, int m)

{
	int i;
	if (e == NULL)
		return;
	simpFree(e->e1, pos, m);
	simpFree(e->e2, pos, m);
	for (i = 0; i < m; i++)
		if (pos[i] == e)
			pos[i] = NULL;
}

static Tree *simplify(Tree *e, Tree **pos, int m)

{
	int i;
	Tree *aux = e;
	switch (e->type)
	{
	case STR: /* nothing to do */
		break;
	case OOR: /* or of two classes is a class, C|eps -> C? */
		e->e1 = simplify(e->e1, pos, m);
		e->e2 = simplify(e->e2, pos, m);
		if ((e->e1->type == STR) && (e->e2->type == STR) &&
			(e->e1->eps == e->e2->eps))
		{
			for (i = 0; i < m; i++)
				if ((pos[i] == e->e1) || (pos[i] == e->e2))
					pos[i] = e;
			e->type = STR;
			e->eps = e->e1->eps;
			free(e->e1);
			e->e1 = NULL;
			free(e->e2);
			e->e2 = NULL;
		}
		else if ((e->e1->type == STR) && e->e1->eps &&
				 (e->e2->type == STR) && !e->e2->eps)
		{
			e->e1->type = QUESTION;
			e->e1->e1 = e->e2;
			e = e->e1;
			free(aux);
		}
		else if ((e->e2->type == STR) && e->e2->eps &&
				 (e->e1->type == STR) && !e->e1->eps)
		{
			e->e2->type = QUESTION;
			e->e2->e1 = e->e1;
			e = e->e2;
			free(aux);
		}
		else if ((e->e2->type == STR) && e->e2->eps &&
				 (e->e1->type == STR) && e->e1->eps)
		{
			freeTree(e->e1);
			e = e->e2;
			free(aux);
		}
		break;
	case CONC: /* conc with epsilon can be removed */
		e->e1 = simplify(e->e1, pos, m);
		e->e2 = simplify(e->e2, pos, m);
		if ((e->e1->type == STR) && e->e1->eps)
		{
			freeTree(e->e1);
			e = e->e2;
			free(aux);
		}
		else if ((e->e2->type == STR) && e->e2->eps)
		{
			freeTree(e->e2);
			e = e->e1;
			free(aux);
		}
		break;
	case STAR: /* ** = *, ?* = *, +* = *, eps* = eps */
		e->e1 = simplify(e->e1, pos, m);
		if ((e->e1->type == PLUS) || (e->e1->type == QUESTION) ||
			(e->e1->type == STAR))
		{
			e = e->e1;
			e->type = STAR;
			e->eps = true;
			free(aux);
		}
		else if ((e->e1->type == STR) && e->e1->eps)
		{
			e = e->e1;
			free(aux);
		}
		break;
	case QUESTION: /* *? = *, ?? = ?, +? = *, eps? = eps */
		e->e1 = simplify(e->e1, pos, m);
		if ((e->e1->type == STAR) || (e->e1->type == PLUS))
		{
			e = e->e1;
			e->type = STAR;
			e->eps = true;
			free(aux);
		}
		else if (e->e1->type == QUESTION)
		{
			e = e->e1;
			free(aux);
		}
		else if ((e->e1->type == STR) && e->e1->eps)
		{
			e = e->e1;
			free(aux);
		}
		break;
	case PLUS: /* remove the + at beginning and ending
			    *+ = *, ?+ = *, ++ = +, eps+ = eps */
		e->e1 = simplify(e->e1, pos, m);
		if ((e->e1->type == STAR) || (e->e1->type == QUESTION))
		{
			e = e->e1;
			e->type = STAR;
			e->eps = true;
			free(aux);
		}
		else if (e->e1->type == PLUS)
		{
			e = e->e1;
			free(aux);
		}
		else if ((e->e1->type == STR) && e->e1->eps)
		{
			e = e->e1;
			free(aux);
		}
		break;
	}
	return e;
}

static Tree *parseLiteral(const char *pat, int *i, Tree **pos)

/* parses a literal pattern */

{
	Tree *e, *r;
	e = (Tree*)malloc(sizeof(Tree));
	e->type = STR;
	e->eps = false;
	e->e1 = e->e2 = NULL;
	e->maskPos = e->firstPos = e->lastPos = NULL;
	pos[*i] = e;
	(*i)++;
	if (!pat[1])
		r = e; /* one letter */
	else
	{
		r = (Tree*)malloc(sizeof(Tree));
		r->type = CONC;
		r->eps = false;
		r->e1 = e;
		r->e2 = parseLiteral(pat + 1, i, pos);
	}
	return r;
}


void setMaskPos(Tree *e, int L)

/* compute maskPos (positions of subexpression) */

{
	switch (e->type)
	{
	case STR:
		e->maskPos = ZERO(createMask(L), L);
		if (!e->eps)
			SET(e->maskPos, e->pos);
		break;
	case STAR:
	case PLUS:
	case QUESTION:
		setMaskPos(e->e1, L);
		e->maskPos = COPY(createMask(L), e->e1->maskPos, L);
		break;
	case OOR:
	case CONC:
		setMaskPos(e->e1, L);
		setMaskPos(e->e2, L);
		e->maskPos = OR(COPY(createMask(L), e->e1->maskPos, L),
						e->e2->maskPos, L);
		break;
	}
}

Tree *parse(const char *pat, int m, Tree **pos, bool simple)

/* captures the structure of the regular expression and
	   maps each character to the place where it has to be put. It then
	   performs the optimizations and returns the real type of the search
	   pattern (simple, extended or regular) and the real number of bits
	   to store it. It returns NULL if there is a syntax error. */

{
	Tree *e;
	int i;
	
	if (m == 0)
		return NULL; /* void pattern */
	for (i = 0; i < m; i++)
		pos[i] = NULL;
	i = 0;
	/* parse the expression */
	if (OptLiteral)
		return parseLiteral(pat, &i, pos);
	e = parseOr(pat, &i, pos);
	if (e == NULL)
		return NULL;
	/* simplify the expression */
	if(simple) e = simplify(e, pos, m);
	return e;
}



std::pair<std::string, std::string> split_pattern(const std::string &pat, int m, Tree** pos, int p_split){
	int i = 0;
	while(i < m){
		if(pos[i] != NULL && pos[i]->pos == p_split){
			while(i < m && (pos[i] == NULL && pat[i] != OpenPar && pat[i] != Or)){
				++i;
			}
			break;
		}
		++i;
	}
	std::pair<std::string, std::string> res;
	res.first = pat.substr(0, i);
	res.second = pat.substr(i);
	return res;

}

