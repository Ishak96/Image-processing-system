#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

typedef struct
{
	int area_label;
	int size;
	int OX, OY;
	float E_x, E_y, E_xy;
	float direction;
	float average_greyscal, avrage_RGB;
	int* H;
} Parameters;

int CEIL = 50;

Parameters init_parameters()
{
	Parameters param = {
		.area_label = -1,
		.size = 0,
		.OX = 0, .OY = 0,
		.E_x = 0, .E_y = 0, .E_xy = 0,
		.direction = 0,
		.average_greyscal = 0, .avrage_RGB = 0,
		.H = ivector(0, 255)
	};

	return param;
}

void free_parameters(Parameters *param)
{
	free_ivector(param->H, 0, 255);
}

// M -> x-axis, N -> y-axis, O -> z-axis
int*** i3D(int M, int N, int O)
{
	int*** A = (int***) malloc(M * sizeof(int**));

	if (A == NULL) {
		fprintf(stderr, "Out of memory");
		exit(0);
	}

	for (int i = 0; i < M; i++) {
		A[i] = (int**) malloc(N * sizeof(int*));

		if (A[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		for (int j = 0; j < N; j++) {
			A[i][j] = (int*)malloc(O * sizeof(int));

	   		if (A[i][j] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
		}
	}

	return A;
}

void free_i3D(int*** A, int M, int N)
{
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++)
			free(A[i][j]);

		free(A[i]);
	}
	free(A);
}

int fileCount(const char* dir) {
	int file_count = 0;
	DIR * dirp;
	struct dirent * entry;

	dirp = opendir(dir); /* There should be error handling after this */	
	int index = 0;
	if(dirp != NULL) {
		while ((entry = readdir(dirp)) != NULL) {
			if(entry->d_type == DT_REG && sscanf(entry->d_name, "lbox%03d.ppm", &index) == 1) {
				file_count++;
			}
		}
	}
	
	closedir(dirp);
	return file_count;	
}

byte** rgb8matrix_to_bmatrix(rgb8** I, int nrl, int nrh, int ncl, int nch)
{
	byte** out = bmatrix(nrl, nrh, ncl, nch);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			out[x][y] = ((int)I[x][y].r + (int)I[x][y].g + (int)I[x][y].b) / 3;
		}
	}

	return out;
}

byte** medianFilter(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq)
{
	int n = strlen(dir_name) + strlen(label) + 3 + strlen(extention);
	const char file_name[n];

	int*** votingBox;

	byte** bOut;

	int som, median;

	for(int i = 1; i <= n_seq ; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir_name, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, nrl, nrh, ncl, nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, *nrl, *nrh, *ncl, *nch);

		if(i == 1) {
			bOut = bmatrix(*nrl, *nrh, *ncl, *nch);
			votingBox = i3D(*nrh, *nch, 256);
		}

		for (int x = *nrl; x < *nrh; x++) {
			for (int y = *ncl; y < *nch; y++) {
				votingBox[x][y][(int) byteI1[x][y]]++;
			}
		}

		free_bmatrix(byteI1, *nrl, *nrh, *ncl, *nch);
		free_rgb8matrix(I1, *nrl, *nrh, *ncl, *nch);
	}

	median = n_seq / 2;

	for (int x = *nrl; x < *nrh; x++) {
		for (int y = *ncl; y < *nch; y++) {
			som = 0;
			for(int z = 0; z < 256; z++) {
				som += votingBox[x][y][z];

				if(som >= median) {
					bOut[x][y] = (byte) z;
					break;
				}
			}
		}
	}

	free_i3D(votingBox, *nrh, *nch);

	return bOut;
}

void minus(byte** I1, byte** I2,long nrl, long nrh, long ncl, long nch)
{

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			int res = (int)I1[x][y] - (int)I2[x][y];
			I2[x][y] = (res >= CEIL) ? (byte) 255 : (byte) 0;
		}
	}	
}

byte** erosion(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
	int same = 0;
	for (int i = 0; i < maskw; i++) {
		for (int j = 0; j < maskh; j++) {
			if(mask[i][j] == 1)
				same++;
		}
	}

   byte** out = bmatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         double acc = 0.0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
        
         out[x][y] = (acc == same) ? (byte) 255 : (byte) 0;
      }
   }

   return out;
}

byte** dilatation(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
   byte** out = bmatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         double acc = 0.0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
        
         out[x][y] = (acc > 0) ? (byte) 255 : (byte) 0;
      }
   }

   return out;
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

Parameters compute_parameters(byte** I, long nrl, long nrh, long ncl, long nch)
{
	Parameters param = init_parameters();


	return param;
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

	for(int i = 1; i <= n_seq; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, nrl, nrh, ncl, nch);

		//traitement
		minus(I, byteI1, nrl, nrh, ncl, nch);
		printf("%d\n", i);
		byte** byteR = cleanImage(byteI1, nrl, nrh, ncl, nch);

		sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
		SavePGM_bmatrix(byteR, nrl, nrh, ncl, nch, file_name_result);
		
		//Free matrix
		free_rgb8matrix(I1, nrl, nrh, ncl, nch);
		
		free_bmatrix(byteI1, nrl, nrh, ncl, nch);
		free_bmatrix(byteR, nrl, nrh, ncl, nch);
	}


	free_bmatrix(I, nrl, nrh, ncl, nch);

	return 0;
}