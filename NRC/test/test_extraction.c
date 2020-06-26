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

extern int CEIL = 40;

int maskH = 3;
int maskW = 3;

int** mask = NULL;

void fill_maskR() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void fill_maskD() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void cleanImage(byte** f, long nrl, long nrh, long ncl, long nch)
{	
	int n = 1;
	
	byte** out = bmatrix(nrl, nrh, ncl, nch);
	
	//debruitage
	fill_maskD();
	n_erosion(f, out, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	n_dilatation(out, f, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	//remplissage
	fill_maskR();

	n_dilatation(f, out, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	n_erosion(out, f, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	free_bmatrix(out, nrl, nrh, ncl, nch);
}

int main(void) {
	
	long nrl, nrh, ncl, nch;
	const char* dir = "../Images/Sequences/Lbox/ppm/";
	const char* label = "lbox";
	const char* extention = "ppm";
	const char* result_dir = "../Images/results/morphoMath/extraxtion/movementDetectiontRefImage/Lbox/";

	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);
	int m = strlen(result_dir) + strlen(label) + 3 + strlen(extention);

	int step = 1;
	int displayNext = step;
	int percent = 0;	

	const char file_name[n];
	const char file_name_result[m];

	int n_seq = fileCount(dir, label);

	byte** I = medianFilter(dir, label, extention, &nrl, &nrh, &ncl, &nch, n_seq);
	SavePGM_bmatrix(I, nrl, nrh, ncl, nch, "../Images/results/morphoMath/extraxtion/lbox_medianFilter.pgm");

	byte** byteI1 = bmatrix(nrl, nrh, ncl, nch);

	mask = imatrix(0, maskH, 0, maskW);

	for(int i = 1; i <= n_seq; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		rgb8matrix_to_bmatrix(I1, byteI1, nrl, nrh, ncl, nch);

		//traitement
		minus(I, byteI1, nrl, nrh, ncl, nch);
		cleanImage(byteI1, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePGM_bmatrix(byteI1, nrl, nrh, ncl, nch, file_name_result);
		
		percent = (100 * (i + 1)) / n_seq;
		if(percent >= displayNext) {
			print_progress(percent, 100);
			displayNext += step;			
		}		

		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
	}

	free_bmatrix(byteI1, nrl, nrh, ncl, nch);
	free_bmatrix(I, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, maskH, 0, maskW);

	return 0;
}