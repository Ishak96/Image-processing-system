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

	int maskE[3][3] = {{1, 1, 1},
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

	//remplissage
	byte **IR = dilatation(ITrou, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		IR = dilatation(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		IR = erosion(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	//debruitage
	n = 8;
	byte **ID = erosion(IBruit, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		ID = erosion(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		ID = dilatation(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}	

	SavePGM_bmatrix(IR, nrl, nrh, ncl, nch, "../Images/results/carre_Trou_remplis.pgm");
	SavePGM_bmatrix(ID, nrl, nrh, ncl, nch, "../Images/results/carre_Bruit_enleve.pgm");
	
	/* Free the byte sum matrix */
	free_bmatrix(ITrou, nrl, nrh, ncl, nch);
	free_bmatrix(IBruit, nrl, nrh, ncl, nch);
	free_bmatrix(IR, nrl, nrh, ncl, nch);
	free_bmatrix(ID, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, 3, 0, 3);

	return 0;
}