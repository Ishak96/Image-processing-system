#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "erosion_dilatation.h"
#include "utils.h"

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

	int n = 1;
	byte **ID = erosion(f, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		ID = erosion(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		ID = dilatation(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}

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
	byte **IR = dilatation(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		IR = dilatation(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		IR = erosion(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	free_imatrix(mask, 0, 3, 0, 3);
	free_bmatrix(ID, nrl, nrh, ncl, nch);

	return IR;
}

int main(void) {
	long nrh, nrl, nch, ncl;

	const char* dir = "../Images/Sequences/Lbox/ppm/";
	const char* result_dir = "../Images/results/morphoMath/movementDetectiontMorphoMath/Lbox/";

	int n = strlen(dir) + 11;
	int m = strlen(result_dir) + 11;

	const char file_name[n];
	const char file_name_result[m];


	int n_seq = fileCount(dir);
	for(int i = 1; i < n_seq-1 ; i++) {
		sprintf(file_name, "%slbox%03d.ppm", dir, i);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, nrl, nrh, ncl, nch);

		//second image
		sprintf(file_name, "%slbox%03d.ppm", dir, i+1);
		rgb8** I2 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		byte** byteI2 = rgb8matrix_to_bmatrix(I2, nrl, nrh, ncl, nch);

		//traitement
		minus(byteI2, byteI1, nrl, nrh, ncl, nch);
		printf("%d\n", i);
		byte** byteR = cleanImage(byteI1, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%slbox%03d.pgm", result_dir, i);
		SavePGM_bmatrix(byteR, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		free_rgb8matrix(I2, nrl, nrh, ncl, nch);
		
		free_bmatrix(byteI1, nrl, nrh, ncl, nch);
		free_bmatrix(byteI2, nrl, nrh, ncl, nch);
		free_bmatrix(byteR, nrl, nrh, ncl, nch);
	}

	return 0;
}