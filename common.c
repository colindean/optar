/* (c) GPL 2007 Karel 'Clock' Kulhavy, Twibright Labs */

#include <stdio.h> /* fprintf */

#include "optar.h"

/* Coordinates don't count with the border - 0,0 is upper left corner of the
 * first cross! */
int is_cross(unsigned x, unsigned y)
{
	x%=CPITCH;
	y%=CPITCH;
	return (x<2*CHALF&&y<2*CHALF);
}

/* Returns the coords relative to the upperloeftmost cross upper left corner
 * pixel! If you have borders, you have to add them! */
void seq2xy(int *x, int *y, unsigned seq)
{
	unsigned rep; /* Repetition - number of narrow strip - wide strip pair,
			 starting with 0 */

	if (seq>=TOTALBITS){
		/* Out of range */
		*x=-1;
		*y=-1;
		return;
	}
	/* We are sure we are in range. Document structure:
	 * - narrow strip (between top row of crosses), height is
	 *   2*CHALF
	 * - wide strip, height is CPITCH-2*CHALF
	 * - the above repeats (YCROSSES-1)-times
	 * - narrow strip 
	 */
	rep=seq/REPPIXELS;
	seq=seq%REPPIXELS;

	*y=REPHEIGHT*rep;
	/* Now seq is sequence in the repetition pair */
	if (seq>=NARROWPIXELS){
		/* Second, wide strip of the pair */
		*y+=NARROWHEIGHT;
		seq-=NARROWPIXELS;
		/* Now seq is sequence in the wide strip */
		*y+=seq/WIDEWIDTH;
		*x=seq%WIDEWIDTH;
	}else{
		/* First, narrow strip of the pair */
		unsigned gap; /* Horizontal gap number */
		*x=2*CHALF;
		*y+=seq/NARROWWIDTH;
		seq%=NARROWWIDTH;
		/* seq is now sequence in the horiz. line */
		gap=seq/GAPWIDTH;
		*x+=gap*CPITCH;
		seq%=GAPWIDTH;
		/* seq is now sequence in the gap */
		*x+=seq;
	}
}

/* Golay codes */
unsigned long golay(unsigned long in)
{
	return golay_codes[in&4095];
}


