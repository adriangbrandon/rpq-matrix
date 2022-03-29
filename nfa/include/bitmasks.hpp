
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

#ifndef BITMASKSINCLUDED
#define BITMASKSINCLUDED

typedef unsigned long mask;

#define ONE ((mask)1)
#define ZEROS ((mask)0)
#define ONES ((mask)~0)

#define W (8*sizeof(mask))
                                                                                
	/* multiword masks. the normal operations are extended, l is
	   the length of the mask in bits */

typedef mask *Mask;
#define maskSize(l) ((l+W-1)/W)

	/* in general op(m,m1) is m <= m op m1; return m */

Mask createMask (int l);
Mask createMasks (int num, int l);
Mask COPY (Mask m, Mask m1, int l);
Mask ZERO (Mask m, int l);
Mask OR (Mask m, Mask m1, int l);
Mask AND (Mask m, Mask m1, int l);
Mask NOT (Mask m, int l);
Mask SET (Mask m, int i);
Mask CLEAR (Mask m, int i);
mask SLICE (Mask m, int from, int len);
boolean EQUAL (Mask m1, Mask m2, int l);
boolean ISSET (Mask m, int i);
boolean ISZERO (Mask m, int l);
int ONEBITS (Mask m, int l);

#endif
