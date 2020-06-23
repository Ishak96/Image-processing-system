#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "erosion_dilatation.h"
#include "extraction.h"

extern int CEIL = 50;

byte** cleanImage(byte** f, long nrl, long nrh, long ncl, long nch)
{
	//debruitage
	float maskD[3][3] = {{0, 1, 0},
						 {1, 0, 1},
						 {0, 1, 0}
						};
	int** mask = imatrix(0, 3, 0, 3);
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask[i][j] = maskD[i][j];
		}
	}

	byte** out11 = bmatrix(nrl, nrh, ncl, nch);
	byte** out12 = bmatrix(nrl, nrh, ncl, nch);

	int n = 1;
	n_erosion(f, out11, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_dilatation(out11, out12, nrl, nrh, ncl, nch, mask, 3, 3, n);


	//remplissage
	float maskE[3][3] = {{0, 1, 0},
						 {1, 1, 1},
						 {0, 1, 0}
						};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask[i][j] = maskE[i][j];
		}
	}

	n = 1;
	n_dilatation(out12, out11, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_erosion(out11, out12, nrl, nrh, ncl, nch, mask, 3, 3, n);

	free_imatrix(mask, 0, 3, 0, 3);
	free_bmatrix(out11, nrl, nrh, ncl, nch);

	return out12;
}

int main(void) {
	
	long nrl, nrh, ncl, nch;
	const char* dir = "../Images/Sequences/Lbox/ppm/";
	const char* label = "lbox";
	const char* extention = "ppm";
	const char* result_dir = "../Images/results/morphoMath/extraxtion/movementDetectiontRefImage/Lbox/";

	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);
	int m = strlen(result_dir) + strlen(label) + 3 + strlen(extention);

	const char file_name[n];
	const char file_name_result[m];

	int n_seq = fileCount(dir);

	byte** I = medianFilter(dir, label, extention, &nrl, &nrh, &ncl, &nch, n_seq);
	SavePGM_bmatrix(I, nrl, nrh, ncl, nch, "../Images/results/morphoMath/extraxtion/lbox_medianFilter.pgm");

	byte** byteI1 = bmatrix(nrl, nrh, ncl, nch);

	for(int i = 1; i <= n_seq; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		rgb8matrix_to_bmatrix(I1, byteI1, nrl, nrh, ncl, nch);

		//traitement
		minus(I, byteI1, nrl, nrh, ncl, nch);
		printf("%d\n", i);
		byte** byteR = cleanImage(byteI1, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePGM_bmatrix(byteR, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		
		free_bmatrix(byteR, nrl, nrh, ncl, nch);
	}

	free_bmatrix(byteI1, nrl, nrh, ncl, nch);
	free_bmatrix(I, nrl, nrh, ncl, nch);

	return 0;
}