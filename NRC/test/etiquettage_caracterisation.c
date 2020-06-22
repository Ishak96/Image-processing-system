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

typedef struct
{
	int area_label;
	int size;
	float OX, OY;
	float E_x, E_y, E_xy;
	float theta;
	float average_greyscal;
	float avrage_R, avrage_G, avrage_B;
	int* H;
} Parameters;

extern int CEIL = 50;

Parameters init_parameters()
{
	Parameters param = {
		.area_label = -1,
		.size = 0,
		.OX = 0, .OY = 0,
		.E_x = 0, .E_y = 0, .E_xy = 0,
		.theta = 0,
		.average_greyscal = 0, 
		.avrage_R = 0, .avrage_G = 0, .avrage_B = 0,
		.H = ivector(0, 255)
	};
	for (int i = 0; i < 255; i++)
		param.H[i] = 0;

	return param;
}

void free_parameters(Parameters *param)
{
	free_ivector(param->H, 0, 255);
}

void dump_parameters(Parameters* param, int size)
{
	printf("########IMAGE#########\n");
	for(int label = 0; label < size; label++) {
		printf("---area label--- : %d\n", param[label].area_label);
		printf("area size : %d\n", param[label].size);
		printf("mean x : %.2f\n", param[label].OX);
		printf("mean y : %.2f\n", param[label].OY);
		printf("mean E_x^2 : %.2f\n", param[label].E_x);
		printf("mean E_y^2 : %.2f\n", param[label].E_y);
		printf("mean E_xy^2 : %.2f\n", param[label].E_xy);
		printf("theta : %.2f\n", param[label].theta);
		printf("average greyscal : %.2f\n", param[label].average_greyscal);
		printf("average RGB [%.2f,%.2f,%.2f] :\n", param[label].avrage_R,
		 										   param[label].avrage_G, param[label].avrage_B);

		printf("[");
		for(int i = 0; i < 254; i++)
			printf("%d, ", param[label].H[i]);

		printf("%d]\n", param[label].H[254]);
	}
	printf("########IMAGE#########\n");
}

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

Parameters* compute_parameters(rgb8** IRGB, int** E, int num_label,
							  long nrl, long nrh, long ncl, long nch)
{
	Parameters* param = (Parameters *) malloc(num_label * sizeof(Parameters));

	for(int label = 1; label <= num_label; label++) {
		Parameters p = init_parameters();
		p.area_label = label;

		ArrayList x_list = NEW_ArrayList();
		ArrayList y_list = NEW_ArrayList();

		for (int y = nrl; y < nrh; y++) {
			for (int x = ncl; x < nch; x++) {
				if(E[y][x] == label) {
					p.size += 1;
					p.OX += x;
					p.OY += y;
					
					int byte = ((int)IRGB[y][x].r + (int)IRGB[y][x].g + (int)IRGB[y][x].b) / 3;
					p.average_greyscal += byte; 

					p.avrage_R += (int)IRGB[y][x].r;
					p.avrage_G += (int)IRGB[y][x].g; 
					p.avrage_B += (int)IRGB[y][x].b;
					p.H[byte] += 1;

					ADD_ArrayList(&y_list, y); ADD_ArrayList(&x_list, x);
				}
			}
		}
		p.OX = p.OX / p.size;
		p.OY = p.OY / p.size;
		p.average_greyscal = p.average_greyscal / p.size;
		p.avrage_R = p.avrage_R / p.size;
		p.avrage_G = p.avrage_G / p.size;
		p.avrage_B = p.avrage_B / p.size;

		for(int i = 0; i < p.size; i++) {
			p.E_x += ((x_list.list[i] - p.OX) * (x_list.list[i] - p.OX)) / p.size;
			p.E_y += ((y_list.list[i] - p.OY) * (y_list.list[i] - p.OY)) / p.size;
			p.E_xy += ((x_list.list[i] - p.OX) * (y_list.list[i] - p.OY)) / p.size;
		}
		p.theta = 0.5 * atan((2 * p.E_x * p.E_y) / (p.E_x - p.E_y));

		param[label-1] = p;
		FREE_ArrayList(&x_list);
		FREE_ArrayList(&y_list);
	}

	return param;
}

void store_parameters(char* file_name, int num_label, Parameters* param)
{

}

void print_progress(size_t count, size_t max)
{
	const char prefix[] = "Progress: [";
	const char suffix[] = "]";
	const size_t prefix_length = sizeof(prefix) - 1;
	const size_t suffix_length = sizeof(suffix) - 1;
	char *buffer = calloc(max + prefix_length + suffix_length + 1, 1);
	size_t i = 0;

	strcpy(buffer, prefix);
	for (; i < max; ++i)
	{
		buffer[prefix_length + i] = i < count ? '#' : ' ';
	}

	strcpy(&buffer[prefix_length + i], suffix);
	printf("\e[1;1H\e[2J");
	printf("\b%c[2K\r%s\n", 27, buffer);
	fflush(stdout);
	free(buffer);
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
		
		dump_parameters(param, nb_label);
		
		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePGM_bmatrix(byteR, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		
		free_bmatrix(byteI1, nrl, nrh, ncl, nch);
		free_bmatrix(byteR, nrl, nrh, ncl, nch);
		free_imatrix(E, nrl, nrh, ncl, nch);
		free(param);

		/*percent = (100 * (i + 1)) / n_seq;
		if(percent >= displayNext) {
			print_progress(percent, 100);
			displayNext += step;			
		}*/		
	}


	free_bmatrix(I, nrl, nrh, ncl, nch);

	return 0;
}