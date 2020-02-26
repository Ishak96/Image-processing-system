#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

int CEIL = 100;

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

byte** temporalAverage(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq)
{
	int n = strlen(dir_name) + strlen(label) + 3 + strlen(extention);
	
	double** out;

	const char file_name[n];

	for(int i = 1; i <= n_seq ; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir_name, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, nrl, nrh, ncl, nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, *nrl, *nrh, *ncl, *nch);

		if(i == 1) {
			out = dmatrix(*nrl, *nrh, *ncl, *nch);			
		}

		for (int x = *nrl; x < *nrh; x++) {
			for (int y = *ncl; y < *nch; y++) {
				out[x][y] = out[x][y] + ((double) byteI1[x][y] / n_seq);
			}
		}

		free_bmatrix(byteI1, *nrl, *nrh, *ncl, *nch);
		free_rgb8matrix(I1, *nrl, *nrh, *ncl, *nch);
	}

	byte** bOut = bmatrix(*nrl, *nrh, *ncl, *nch);			

	for (int x = *nrl; x < *nrh; x++) {
		for (int y = *ncl; y < *nch; y++) {
			bOut[x][y] = (byte) abs(out[x][y]);
		}
	}

	return bOut;
}

void minus(byte** I1, byte** I2,long nrl, long nrh, long ncl, long nch)
{

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			int res = abs((int)I1[x][y] - (int)I2[x][y]);
			I1[x][y] = (res > CEIL) ? (byte) 255 : (byte) 0;
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
	float maskE[1][1] = {{1}};
   	float** mask = imatrix(0, 1, 0, 1);

	for (int i = 0; i < 1; i++) {
		for (int j = 0; j < 1; j++) {
			mask[i][j] = maskE[i][j];
		}
	}

	int n = 5;

	//remplissage
	byte **IR = dilatation(f, nrl, nrh, ncl, nch, mask, 1, 1);
	for(int i = 1; i < n; i++) {
		IR = dilatation(IR, nrl, nrh, ncl, nch, mask, 1, 1);
	}

	for(int i = 0; i < n; i++) {
		IR = erosion(IR, nrl, nrh, ncl, nch, mask, 1, 1);
	}

	free_imatrix(mask, 0, 1, 0, 1);

	//debruitage
	float maskD[2][2] = {{1, 1},
						 {1, 1}
						};
   	mask = imatrix(0, 2, 0, 2);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			mask[i][j] = maskD[i][j];
		}
	}

	n = 12;
	byte **ID = erosion(IR, nrl, nrh, ncl, nch, mask, 2, 2);
	for(int i = 1; i < n; i++) {
		ID = erosion(ID, nrl, nrh, ncl, nch, mask, 2, 2);
	}

	for(int i = 0; i < n; i++) {
		ID = dilatation(ID, nrl, nrh, ncl, nch, mask, 2, 2);
	}

	free_imatrix(mask, 0, 2, 0, 2);
	//free_imatrix(IR, nrl, nrh, ncl, nch);

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

	byte** I = temporalAverage(dir, label, extention, &nrl, &nrh, &ncl, &nch, n_seq);
	
	for(int i = 1; i <= n_seq; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);
		byte** byteI1 = rgb8matrix_to_bmatrix(I1, nrl, nrh, ncl, nch);

		//traitement
		minus(byteI1, I, nrl, nrh, ncl, nch);
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