#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"

#define LAMBDA 0.00005

rgb8** gradient(byte** I, long nrl, long nrh, long ncl, long nch)
{
	rgb8** harris_i = rgb8matrix(nrl, nrh, ncl, nch);

	float maskd1[3][3] = {{-1, 0, 1},
      					 { -2, 0, 2},
      					 { -1, 0, 1}
   						 };

	float maskd2[3][3] = {{-1, -2, -1},
      					 {  0,  0,  0},
      					 {  1,  2,  1}
   						 };

	float maskg[3][3] = {{1, 1, 1},
      					 {1, 0, 1},
      					 {1, 1, 1}
   						};

	float** mask1 = imatrix(0, 3, 0, 3);
	float** mask2 = imatrix(0, 3, 0, 3);
	float** mask3 = imatrix(0, 3, 0, 3);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask1[i][j] = maskd1[i][j];
			mask2[i][j] = maskd2[i][j];
			mask3[i][j] = maskg[i][j];
		}
	}

	/* Calculate convolutions for the two masks */
	int** ix = conv2(I, nrl, nrh, ncl, nch, mask1, 3, 3); // horizontal
	int** iy = conv2(I, nrl, nrh, ncl, nch, mask2, 3, 3); // vertical

	int** ixy = multiply(ix, iy, nrl, nrh, ncl, nch);

	int** ixy_g = gaussian_filter(ixy, nrl, nrh, ncl, nch, mask3, 3, 3);

	apply(ix, nrl, nrh, ncl, nch, pow2);
	apply(iy, nrl, nrh, ncl, nch, pow2);

	int** ixx = gaussian_filter(ix, nrl, nrh, ncl, nch, mask3, 3, 3);
	int** iyy = gaussian_filter(iy, nrl, nrh, ncl, nch, mask3, 3, 3);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			harris_i[x][y].r = I[x][y]; harris_i[x][y].g = I[x][y]; harris_i[x][y].b = I[x][y];
			
			float a = ix[x][y]*iyy[x][y] + iy[x][y]*ixx[x][y] - 2*ixy[x][y]*ixy_g[x][y];
			float b = ixx[x][y] + iyy[x][y];

			float res = b == 0 ? 0 : a / b; 

			if(res > 3500) {
				harris_i[x][y].r = 255;
				harris_i[x][y].g = 0; 
				harris_i[x][y].b = 0;
			}
		}
	}

	free_imatrix(ix, nrl, nrh, ncl, nch);
	free_imatrix(iy, nrl, nrh, ncl, nch);
	free_imatrix(ixx, nrl, nrh, ncl, nch);
	free_imatrix(iyy, nrl, nrh, ncl, nch);
	free_imatrix(ixy, nrl, nrh, ncl, nch);
	free_imatrix(ixy_g, nrl, nrh, ncl, nch);

	free_imatrix(mask1, 0, 3, 0, 3);
	free_imatrix(mask2, 0, 3, 0, 3);
	free_imatrix(mask3, 0, 3, 0, 3);

	return harris_i;
}

rgb8** harris(byte** I, long nrl, long nrh, long ncl, long nch)
{
	rgb8** harris_i = rgb8matrix(nrl, nrh, ncl, nch);

	float maskd1[3][3] = {{-1, 0, 1},
      					 { -2, 0, 2},
      					 { -1, 0, 1}
   						 };

	float maskd2[3][3] = {{-1, -2, -1},
      					 {  0,  0,  0},
      					 {  1,  2,  1}
   						 };

	float maskg[3][3] = {{1, 2, 1},
      					 {2, 4, 2},
      					 {1, 2, 1}
   						};

	float** mask1 = imatrix(0, 3, 0, 3);
	float** mask2 = imatrix(0, 3, 0, 3);
	float** mask3 = imatrix(0, 3, 0, 3);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask1[i][j] = maskd1[i][j];
			mask2[i][j] = maskd2[i][j];
			mask3[i][j] = maskg[i][j];
		}
	}

	/* Calculate convolutions for the two masks */
	int** ix = conv2(I, nrl, nrh, ncl, nch, mask1, 3, 3); // horizontal
	int** iy = conv2(I, nrl, nrh, ncl, nch, mask2, 3, 3); // vertical
	
	int** ixy = multiply(ix, iy, nrl, nrh, ncl, nch);
	ixy = gaussian_filter(ixy, nrl, nrh, ncl, nch, mask3, 3, 3);
	
	apply(ix, nrl, nrh, ncl, nch, pow2);
	apply(iy, nrl, nrh, ncl, nch, pow2);
	apply(ixy, nrl, nrh, ncl, nch, pow2);

	int** ixx = gaussian_filter(ix, nrl, nrh, ncl, nch, mask3, 3, 3);
	int** iyy = gaussian_filter(iy, nrl, nrh, ncl, nch, mask3, 3, 3);

	int** ixy_2 = multiply(ixx, iyy, nrl, nrh, ncl, nch);

	int** sum_ixiy_2 = sum(ixx, iyy, nrl, nrh, ncl, nch);
	apply(sum_ixiy_2, nrl, nrh, ncl, nch, pow2);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			harris_i[x][y].r = I[x][y]; harris_i[x][y].g = I[x][y]; harris_i[x][y].b = I[x][y];
			
			float res = ixy_2[x][y] - ixy[x][y] - LAMBDA * sum_ixiy_2[x][y];
			
			if(res > 3500) {
				harris_i[x][y].r = 255;
				harris_i[x][y].g = 0; 
				harris_i[x][y].b = 0;
			}
		}
	}

	free_imatrix(ix, nrl, nrh, ncl, nch);
	free_imatrix(iy, nrl, nrh, ncl, nch);
	free_imatrix(ixx, nrl, nrh, ncl, nch);
	free_imatrix(iyy, nrl, nrh, ncl, nch);
	free_imatrix(ixy, nrl, nrh, ncl, nch);
	free_imatrix(ixy_2, nrl, nrh, ncl, nch);
	free_imatrix(sum_ixiy_2, nrl, nrh, ncl, nch);

	free_imatrix(mask1, 0, 3, 0, 3);
	free_imatrix(mask2, 0, 3, 0, 3);
	free_imatrix(mask3, 0, 3, 0, 3);

	return harris_i;
}

int main(void) {
	long nrh, nrl,
         nch, ncl;
	byte **I;

	I = LoadPGM_bmatrix("../Images/src/cubesx3.pgm", &nrl, &nrh, &ncl, &nch);

	rgb8** harris_i = harris(I, nrl, nrh, ncl, nch);

	SavePPM_rgb8matrix(harris_i, nrl, nrh, ncl, nch, "../Images/results/corner_detection/gradient_corner_detection_cubesx3.ppm");

	/* Free the byte sum matrix */
	free_bmatrix(I, nrl, nrh, ncl, nch);
	free_rgb8matrix(harris_i, nrl, nrh, ncl, nch);
	
	return 0;
}