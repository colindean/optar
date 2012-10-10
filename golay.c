/* (c) GPL 2007 Karel 'Clock' Kulhavy, Twibright Labs */

#include <stdio.h>
#include <assert.h>

#include "parity.h"

int dodecahedron[12][5]={
	/* For each dodecahedron face (number in the comment, 1-12) there
	 * is a list of the 5 adjacent faces (1-12). See golay.svg for
	 * a drawing. */
	{/*  1 */   2,  3,  4,  5,  6},
	{/*  2 */   1,  3,  6,  7,  8},
	{/*  3 */   1,  2,  4,  8,  9},
	{/*  4 */   1,  3,  5,  9, 10},
	{/*  5 */   1,  4,  6, 10, 11},
	{/*  6 */   1,  2,  5,  7, 11},
	{/*  7 */   2,  6,  8, 11, 12},
	{/*  8 */   2,  3,  7,  9, 12},
	{/*  9 */   3,  4,  8, 10, 12},
	{/* 10 */   4,  5,  9, 11, 12},
	{/* 11 */   5,  6,  7, 10, 12},
	{/* 12 */   7,  8,  9, 10, 11}
};

unsigned parities[12];

int main(int argc, char ** argv)
{
	unsigned mask, p, f;
	unsigned input;
	unsigned prty; /* parity */

	for (p=0;p<12;p++){
		mask=0xfff; /* All dodecahedron faces */
		for (f=0;f<5;f++)
			mask^=1U<<(dodecahedron[p][f]-1);
		parities[p]=mask;
	}

	printf("unsigned long golay_codes[4096]={\n");

	for (input=0;input<4096;input++){
		unsigned n_ones;
		unsigned long codeword;

		prty=0;
		for (p=0;p<12;p++){
			prty<<=1;
			prty|=parity(input&parities[p]);
		}
		codeword=((unsigned long)input<<12)|prty;
		n_ones=ones(codeword);
		printf((input==4095?"0x%06lx\n":"0x%06lx,\n"),codeword);
		assert(n_ones==0
			||n_ones==8
			||n_ones==12
			||n_ones==16
			||n_ones==24);
	}

	printf("};\n");
	return 0;
}
