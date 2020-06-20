#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"

extern int CEIL = 60;

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

	//debruitage
	n = 12;
	byte **ID = erosion(IR, nrl, nrh, ncl, nch, mask, 1, 1);
	for(int i = 1; i < n; i++) {
		ID = erosion(ID, nrl, nrh, ncl, nch, mask, 1, 1);
	}

	for(int i = 0; i < n; i++) {
		ID = dilatation(ID, nrl, nrh, ncl, nch, mask, 1, 1);
	}

	free_imatrix(mask, 0, 1, 0, 1);
	free_imatrix(IR, nrl, nrh, ncl, nch);

	return ID;
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
		minus(byteI1, byteI2, nrl, nrh, ncl, nch);
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