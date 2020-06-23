#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "erosion_dilatation.h"
#include "utils.h"

int main(void) {
	long nrh, nrl, nch, ncl;

	int maskE[3][3] = { {1, 1, 1},
      					{1, 1, 1},
      					{1, 1, 1}
   					  };
   int** mask = imatrix(0, 3, 0, 3);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask[i][j] = maskE[i][j];
		}
	}
	
	int n = 6;

	byte** ITrou = LoadPGM_bmatrix("../Images/Test/carreTrou.pgm", &nrl, &nrh, &ncl, &nch);
	byte** IBruit = LoadPGM_bmatrix("../Images/Test/carreBruit.pgm", &nrl, &nrh, &ncl, &nch);

	byte** out11 = bmatrix(nrl, nrh, ncl, nch);
	byte** out12 = bmatrix(nrl, nrh, ncl, nch);

	byte** out21 = bmatrix(nrl, nrh, ncl, nch);
	byte** out22 = bmatrix(nrl, nrh, ncl, nch);

	//remplissage
	n_dilatation(ITrou, out11, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_erosion(out11, out12, nrl, nrh, ncl, nch, mask, 3, 3, n);

	//debruitage
	n = 8;
	n_erosion(IBruit, out21, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_dilatation(out21, out22, nrl, nrh, ncl, nch, mask, 3, 3, n);

	SavePGM_bmatrix(out12, nrl, nrh, ncl, nch, "../Images/results/carre_Trou_remplis.pgm");
	SavePGM_bmatrix(out22, nrl, nrh, ncl, nch, "../Images/results/carre_Bruit_enleve.pgm");
	
	/* Free the byte sum matrix */
	free_bmatrix(ITrou, nrl, nrh, ncl, nch);
	free_bmatrix(IBruit, nrl, nrh, ncl, nch);
	free_bmatrix(out11, nrl, nrh, ncl, nch);
	free_bmatrix(out12, nrl, nrh, ncl, nch);
	free_bmatrix(out21, nrl, nrh, ncl, nch);
	free_bmatrix(out22, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, 3, 0, 3);

	return 0;
}