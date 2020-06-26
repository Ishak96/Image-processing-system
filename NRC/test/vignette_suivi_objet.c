#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>
#include <gtk/gtk.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "erosion_dilatation.h"
#include "extraction.h"
#include "labelling.h"
#include "etiquettage_caracterisation.h"
#include "track.h"

extern int CEIL = 25;

int** mask = NULL;

int maskH = 2;
int maskW = 2;

void fill_maskR() {
	mask[0][0] = 1; mask[0][1] = 1; 
	mask[1][0] = 1; mask[1][1] = 1;
}

void fill_maskD() {
	mask[0][0] = 1; mask[0][1] = 1; 
	mask[1][0] = 0; mask[1][1] = 1;
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
	const char* ref_image_file = "../Images/results/morphoMath/extraxtion/lbox_medianFilter.pgm";

	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);
	int m = strlen(result_dir) + strlen(label) + 3 + strlen(extention);

	const char file_name[n];
	const char file_name_result[m];

	int n_seq = fileCount(dir, label);
	
	int step = 1;
	int displayNext = step;
	int percent = 0;	

	byte** I;

	Vignette vignette;
	
	int vignette_size = 30;

	if(access(ref_image_file, F_OK) != -1) {
		I = LoadPGM_bmatrix(ref_image_file, &nrl, &nrh, &ncl, &nch);
	} else {
		I = medianFilter(dir, label, extention, &nrl, &nrh, &ncl, &nch, n_seq);
		SavePGM_bmatrix(I, nrl, nrh, ncl, nch, ref_image_file);
	}

	byte** byteI1 = bmatrix(nrl, nrh, ncl, nch);
	mask = imatrix(0, maskH, 0, maskW);
	int** harris_i = imatrix(nrl, nrh, ncl, nch);
	int** E = imatrix(nrl, nrh, ncl, nch);
	rgb8** I1;

	int track_label = 1;

	for(int i = 1; i <= n_seq; i++) {

		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		
		int nb_label;

		//traitement
		rgb8_minus(I1, I, byteI1, nrl, nrh, ncl, nch);
		
		cleanImage(byteI1, nrl, nrh, ncl, nch);

		LOOKUP_TABLE_LABELLING(byteI1, E, &nb_label, nrl, nrh, ncl, nch);

		if(i == 1){
			harris(harris_i, byteI1, nrl, nrh, ncl, nch);
			init_system(&vignette, I1, harris_i, E, 
						track_label, vignette_size, nrl, nrh, ncl, nch);
			dump_vignette(vignette);
		}

		track_label = track(E, I1, &vignette, nrl, nrh, ncl, nch);
		if(track_label != 0)
			extract_zone(I1, byteI1, vignette, 0, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePPM_rgb8matrix(I1, nrl, nrh, ncl, nch, file_name_result);

		free_rgb8matrix(I1, nrl, nrh, ncl, nch);

		percent = (100 * (i + 1)) / n_seq;
		if(percent >= displayNext) {
			print_progress(percent, 100);
			displayNext += step;			
		}
	}

	//Free matrix
	free_bmatrix(byteI1, nrl, nrh, ncl, nch);
	free_imatrix(harris_i, nrl, nrh, ncl, nch);
	free_imatrix(E, nrl, nrh, ncl, nch);
	destroy_mask();
	free_imatrix(mask, 0, maskH, 0, maskW);
	free_bmatrix(I, nrl, nrh, ncl, nch);

	return 0;
}