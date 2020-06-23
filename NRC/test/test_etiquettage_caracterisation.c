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
#include "labelling.h"
#include "etiquettage_caracterisation.h"

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
	//remplissage
	byte **IR = dilatation(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		IR = dilatation(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		IR = erosion(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	free_imatrix(mask, 0, 3, 0, 3);
	free_imatrix(ID, nrl, nrh, ncl, nch);

	return IR;
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
	
	int step = 1;
	int displayNext = step;
	int percent = 0;	

	byte** I = medianFilter(dir, label, extention, &nrl, &nrh, &ncl, &nch, n_seq);
	SavePGM_bmatrix(I, nrl, nrh, ncl, nch, "../Images/results/morphoMath/extraxtion/lbox_medianFilter.pgm");

	for(int i = 1; i <= n_seq; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, nrl, nrh, ncl, nch);
		int** E;
		int nb_label;

		//traitement
		minus(I, byteI1, nrl, nrh, ncl, nch);
		
		byte** byteR = cleanImage(byteI1, nrl, nrh, ncl, nch);

		E = LOOKUP_TABLE_LABELLING(byteR, &nb_label, nrl, nrh, ncl, nch);

		Parameters* param = compute_parameters(I1, E, nb_label,
							  				   nrl, nrh, ncl, nch);
		
		store_parameters(i, nb_label, param);
		
		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePGM_bmatrix(byteR, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		
		free_bmatrix(byteI1, nrl, nrh, ncl, nch);
		free_bmatrix(byteR, nrl, nrh, ncl, nch);
		free_imatrix(E, nrl, nrh, ncl, nch);
		free(param);

		percent = (100 * (i + 1)) / n_seq;
		if(percent >= displayNext) {
			print_progress(percent, 100);
			displayNext += step;			
		}
	}


	free_bmatrix(I, nrl, nrh, ncl, nch);

	return 0;
}