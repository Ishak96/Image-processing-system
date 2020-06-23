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

int** mask = NULL;

void fill_maskR() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void fill_maskD() {
	mask[0][0] = 0; mask[0][1] = 1; mask[0][2] = 0;
	mask[1][0] = 0; mask[1][1] = 1; mask[1][2] = 0;
	mask[2][0] = 0; mask[2][1] = 1; mask[2][2] = 0;
}

void cleanImage(byte** f, long nrl, long nrh, long ncl, long nch)
{	
	int n = 1;
	
	byte** out = bmatrix(nrl, nrh, ncl, nch);
	
	//remplissage
	fill_maskR();

	n_dilatation(f, out, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_erosion(out, f, nrl, nrh, ncl, nch, mask, 3, 3, n);	
	
	//debruitage
	n = 2;
	fill_maskD();
	n_erosion(f, out, nrl, nrh, ncl, nch, mask, 3, 3, n);

	n_dilatation(out, f, nrl, nrh, ncl, nch, mask, 3, 3, n);

	free_bmatrix(out, nrl, nrh, ncl, nch);
}

int main(void) {
	long nrh, nrl, nch, ncl;

	const char* dir = "../Images/Sequences/Lbox/ppm/";
	const char* result_dir = "../Images/results/morphoMath/movementDetectiontMorphoMath/Lbox/";

	int n = strlen(dir) + 11;
	int m = strlen(result_dir) + 11;

	const char file_name[n];
	const char file_name_result[m];

	int step = 1;
	int displayNext = step;
	int percent = 0;	

	byte** byteI1 = NULL;
	byte** byteI2 = NULL;
	mask = imatrix(0, 3, 0, 3);

	int n_seq = fileCount(dir);
	for(int i = 1; i < n_seq-1 ; i++) {
		sprintf(file_name, "%slbox%03d.ppm", dir, i);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);

		if(i == 1) {
			byteI1 = bmatrix(nrl, nrh, ncl, nch);
			byteI2 = bmatrix(nrl, nrh, ncl, nch);
		}

		rgb8matrix_to_bmatrix(I1, byteI1, nrl, nrh, ncl, nch);

		//second image
		sprintf(file_name, "%slbox%03d.ppm", dir, i+1);
		rgb8** I2 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);

		rgb8matrix_to_bmatrix(I2, byteI2, nrl, nrh, ncl, nch);

		//traitement
		minus(byteI2, byteI1, nrl, nrh, ncl, nch);

		//cleanImage(byteI1, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%slbox%03d.pgm", result_dir, i);
		SavePGM_bmatrix(byteI1, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		free_rgb8matrix(I2, nrl, nrh, ncl, nch);

		percent = (100 * (i + 1)) / n_seq;
		if(percent >= displayNext) {
			print_progress(percent, 100);
			displayNext += step;			
		}		
	}

	free_bmatrix(byteI1, nrl, nrh, ncl, nch);
	free_bmatrix(byteI2, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, 3, 0, 3);

	return 0;
}