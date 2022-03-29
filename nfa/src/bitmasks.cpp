
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

#include "../include/basics.hpp"

	/* multiword masks */

Mask createMask (int l)

   { return (Mask) malloc (maskSize(l) * sizeof(mask));
   }

Mask createMasks (int num, int l)

   { return (Mask) malloc (num * maskSize(l) * sizeof(mask));
   }

Mask COPY (Mask m, Mask m1, int l)

  { int L = maskSize(l);
    int j;
    for (j=0;j<L; j++) m[j] = m1[j];
    return m;
  }

Mask ZERO (Mask m, int l)

  { int L = maskSize(l);
    int j;
    for (j=0;j<L; j++) m[j] = ZEROS;
    return m;
  }

Mask OR (Mask m, Mask m1, int l)

  { int L = maskSize(l);
    int j;
    for (j=0;j<L; j++) m[j] |= m1[j];
    return m;
  }

Mask AND (Mask m, Mask m1, int l)

  { int L = maskSize(l);
    int j;
    for (j=0;j<L; j++) m[j] &= m1[j];
    return m;
  }

Mask NOT (Mask m, int l)

  { int L = maskSize(l);
    int j;
    for (j=0;j<L; j++) m[j] = ~m[j];
    return m;
  }

boolean EQUAL (Mask m1, Mask m2, int l)

  { int L = maskSize(l);
    int j;
    for (j=0; j<L; j++)
       if (m1[j] != m2[j]) return false;
    return true;
  }

boolean ISSET (Mask m, int i)

  { return (m[i/W] & (ONE << (i % W))) != 0;
  }

boolean ISZERO (Mask m, int l)

  { int L = maskSize(l);
    int j;
    for (j=0; j<L; j++) if (m[j]) return false;
    return true;
  }

int ONEBITS (Mask m, int l)

  { int j,n=0;
    for (j=0; j<l; j++) 
	if (ISSET(m,j)) n++;
    return n;
  }

Mask SET (Mask m, int i)

  { m[i/W] |= ONE << (i % W);
    return m;
  }

Mask CLEAR (Mask m, int i)

  { m[i/W] &= ~ (ONE << (i % W));
    return m;
  }

mask SLICE (Mask m, int from, int len)

  { int i = from / W;
    mask msk = (len == W) ? ONES : ((ONE<<len)-ONE);
    from %= W;
    return (m[i] >> from) & msk;
  }

