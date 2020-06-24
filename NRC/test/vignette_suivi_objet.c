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

typedef struct {
	
	int OX, OY;
	int min_x, max_x; 
	int min_y, max_y;
} Pose;

typedef struct {
	
	int size;
	Pose pose;
	float average_r, average_g, average_b;
	int area_size;
} Vignette;

extern int CEIL = 25;

int** mask = NULL;

int maskH = 3;
int maskW = 3;

void fill_maskR() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void fill_maskD() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 0;
	mask[1][0] = 0; mask[1][1] = 1; mask[1][2] = 0;
	mask[2][0] = 0; mask[2][1] = 1; mask[2][2] = 1;
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

void rgb8_minus(rgb8** I, byte** Iref, byte** I2, long nrl, long nrh, long ncl, long nch)
{
	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			int value = ((int)I[x][y].r + (int)I[x][y].g + (int)I[x][y].b) / 3;
			int res = (int)Iref[x][y] - value;
			I2[x][y] = (res >= CEIL) ? (byte) 255 : (byte) 0;
		}
	}	
}

void init_system(Vignette* vignette, rgb8** I, int** harris_i, int** E, 
				 int track_label, int vignette_size, long nrl, long nrh, long ncl, long nch)
{
	vignette->pose.OY = -INFINITY;	

	int N = 0;
	int sum_r = 0, sum_g = 0, sum_b = 0;

	for (int y = nrl; y < nrh; y++) {
		for (int x = ncl; x < nch; x++) {
			if (E[y][x] == track_label) {
            	sum_r += (int)I[y][x].r;
            	sum_g += (int)I[y][x].g;
            	sum_b += (int)I[y][x].b;
            	N += 1;
            
				if(y > vignette->pose.OY && harris_i[y][x]) {
					vignette->pose.OX = x;
					vignette->pose.OY = y;
				}
			}
		}
	}

	vignette->average_r = (float)sum_r / (float)N;
	vignette->average_g = (float)sum_g / (float)N;
	vignette->average_b = (float)sum_b / (float)N;
	vignette->size = vignette_size;
	vignette->area_size = N;
	vignette->pose.min_x = vignette->pose.max_x = 0;
	vignette->pose.min_y = vignette->pose.max_y = 0;
}

void compute_param(rgb8** I, int** E, int label, float* next_r, float* next_g, float* next_b,
				   int* next_s, int* min_x, int* max_x, int* min_y, int* max_y,
				   long nrl, long nrh, long ncl, long nch) {
	int N = 0;
	int sum_r = 0, sum_g = 0, sum_b = 0;

	for (int y = nrl; y < nrh; y++) {
		for (int x = ncl; x < nch; x++) {
			if (E[y][x] == label) {
				
				if(x > *max_x) *max_x = x;
				if(x < *min_x) *min_x = x;
            
				if(y > *max_y) *max_y = y;
				if(y < *min_y) *min_y = y;

            	sum_r += (int)I[y][x].r;
            	sum_g += (int)I[y][x].g;
            	sum_b += (int)I[y][x].b;
            	N += 1;
            }
        }
    }

    *next_s = N;
	*next_r = (float)sum_r / (float)N;
	*next_g = (float)sum_g / (float)N;
	*next_b = (float)sum_b / (float)N;
}

void dump_vignette(Vignette vignette)
{
	printf("size : %d\n", vignette.size);
	printf("average rgb [%.2f,%.2f,%.2f]\n", vignette.average_r
										   , vignette.average_g
										   , vignette.average_b);
	printf("position : (%d %d)\n", vignette.pose.OY, vignette.pose.OX);	
}

int track(int** E, rgb8** I1, Vignette* vignette, long nrl, long nrh, long ncl, long nch)
{
	int new_label = 0;
	int break_c = 0;
	int err = 1000;

	float new_r = 0, new_g = 0, new_b = 0;
	int new_x = vignette->pose.OX, new_y = vignette->pose.OY;
	
	
	//the current avrage greyscal
	float curr_r = vignette->average_r;
	float curr_g = vignette->average_g;
	float curr_b = vignette->average_b;
	int curr_s = vignette->area_size;
	
	int curr_min_x = vignette->pose.min_x;
	int curr_max_x = vignette->pose.max_x;
	int curr_min_y = vignette->pose.min_y;
	int curr_max_y = vignette->pose.max_y;

	for (int u = 0; u < vignette->size && !break_c; u++) {
		for (int v = 0; v < vignette->size && !break_c; v++) {
				
			// The thumbnail centred on the point of interest
			int ny = vignette->pose.OY + u - ((int) (vignette->size / 2));
			int nx = vignette->pose.OX + v - ((int) (vignette->size / 2));

			if(ny >= nrl && ny < nrh && nx >= ncl && nx < nch) {
					
				int label = E[ny][nx];
				if(label != 0) {
					float next_r = 0, next_g = 0, next_b = 0;
					int next_s = 0;
					int next_min_x = INFINITY, next_max_x = -INFINITY, 
						next_min_y = INFINITY, next_max_y = -INFINITY;

					compute_param(I1, E, label, &next_r, &next_g, &next_b,
								  &next_s, &next_min_x, &next_max_x, &next_min_y, &next_max_y,
								  nrl, nrh, ncl, nch);

					int new_err =  (int) (100 * ((abs(next_r-curr_r) / curr_r
									+ abs(next_g-curr_g) / curr_g
									+ abs(next_b-curr_b) / curr_b
									+ abs(next_s-curr_s) / curr_s) / 4));

					if(new_err < err) {
						err = new_err;
						new_label = label;
						
						new_y = ny;
						new_x = nx;
						
						curr_min_x = next_min_x;
						curr_max_x = next_max_x;
						curr_min_y = next_min_y;
						curr_max_y = next_max_y;

						if(err == 0)
							break_c = 1;
					}
				}
			}
		}
	}
	
	vignette->pose.OY = new_y;
	vignette->pose.OX = new_x;

	vignette->pose.min_x = curr_min_x;
	vignette->pose.max_x = curr_max_x;
	vignette->pose.min_y = curr_min_y;
	vignette->pose.max_y = curr_max_y;	

	return new_label;
}

void extract_zone(rgb8** I, byte** byteI1, Vignette vignette,
				  int rgb, long nrl, long nrh, long ncl, long nch)
{
	int min_x = vignette.pose.min_x, min_y = vignette.pose.min_y; 
	int max_x = vignette.pose.max_x, max_y = vignette.pose.max_y;
	
	if(!rgb) {
		for (int y = nrl; y < nrh; y++) {
			for (int x = ncl; x < nch; x++) {
				I[y][x].r = I[y][x].g = I[y][x].b = (int) byteI1[y][x];
			}
		}
	}
	
	for (int y = min_y; y < max_y; y++) {
		for (int x = min_x; x < max_x; x++) {
			if(x >= max_x - 1 || y >= max_y - 1 || x <= min_x  || y <= min_y) {
				I[y][x].r = 255; I[y][x].g = 0; I[y][x].b = 0;
			}
		}
	}
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

	int n_seq = fileCount(dir);
	
	int step = 1;
	int displayNext = step;
	int percent = 0;	

	byte** I;

	Vignette vignette;
	
	int vignette_size = 25;

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
			extract_zone(I1, byteI1, vignette, 1, nrl, nrh, ncl, nch);

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