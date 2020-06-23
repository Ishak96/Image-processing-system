#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "corner_detection.h"
#include "utils.h"

float** mask1 = NULL;
float** mask2 = NULL;
float** mask3 = NULL;

void create_mask(int harris_gradient) 
{
	if(mask1 == NULL) {
		mask1 = imatrix(0, 3, 0, 3);
		mask1[0][0] = -1; mask1[0][1] = 0; mask1[0][2] = 1;
		mask1[1][0] = -2; mask1[1][1] = 0; mask1[1][2] = 2;
		mask1[2][0] = -1; mask1[2][1] = 0; mask1[2][2] = 1;
	}

	if(mask2 == NULL) {
		mask2 = imatrix(0, 3, 0, 3);
		mask2[0][0] = -1; mask2[0][1] = -2; mask2[0][2] = -1;
		mask2[1][0] = 0; mask2[1][1] = 0; mask2[1][2] = 0;
		mask2[2][0] = 1; mask2[2][1] = 2; mask2[2][2] = 1;
	}

	if(mask3 == NULL) {
		mask3 = imatrix(0, 3, 0, 3);
		mask3[0][0] = 1;  mask3[0][2] = 1;
		if(harris_gradient){
			mask3[0][1] = 2;
			mask3[1][0] = 2; mask3[1][1] = 4; mask3[1][2] = 2;
			mask3[2][1] = 2;
		} else {
			mask3[0][1] = 1;
			mask3[1][0] = 1; mask3[1][1] = 0; mask3[1][2] = 1;
			mask3[2][1] = 1;			
		}
		mask3[2][0] = 1; mask3[2][2] = 1;
	}	
}

void destroy_mask()
{
	if(mask1 != NULL && mask2 != NULL && mask3 != NULL) {
		free_imatrix(mask1, 0, 3, 0, 3);
		free_imatrix(mask2, 0, 3, 0, 3);
		free_imatrix(mask3, 0, 3, 0, 3);
	}
}

int** gradient(byte** I, long nrl, long nrh, long ncl, long nch)
{
	create_mask(0);

	int** grad_i = imatrix(nrl, nrh, ncl, nch);

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

			float a = ix[x][y]*iyy[x][y] + iy[x][y]*ixx[x][y] - 2*ixy[x][y]*ixy_g[x][y];
			float b = ixx[x][y] + iyy[x][y];

			float res = b == 0 ? 0 : a / b; 

			if(res > 3500) 
				grad_i[x][y] = 1;
		}
	}

	free_imatrix(ix, nrl, nrh, ncl, nch);
	free_imatrix(iy, nrl, nrh, ncl, nch);
	free_imatrix(ixx, nrl, nrh, ncl, nch);
	free_imatrix(iyy, nrl, nrh, ncl, nch);
	free_imatrix(ixy, nrl, nrh, ncl, nch);
	free_imatrix(ixy_g, nrl, nrh, ncl, nch);

	return grad_i;
}

void harris(int** harris_i, byte** I, long nrl, long nrh, long ncl, long nch)
{
	create_mask(1);
	
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

	int local_x = 0;
	int local_y = 0;
	float local_max = -10000;
	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			
			float res = (ixx[x][y] * iyy[x][y]) - ixy[x][y] - LAMBDA * (ixx[x][y] + iyy[x][y]) * (ixx[x][y] + iyy[x][y]);
			harris_i[x][y] = 0;
			
			if(res < 0) {
				if(local_max > 0) {
					harris_i[local_x][local_y] = 1;
					local_max = -10000;
				}
			} else {
				if(res > local_max && res > 3500) {
					local_max = res;
					local_x = x;
					local_y = y;
				}
			}
		}
	}

	free_imatrix(ixx, nrl, nrh, ncl, nch);
	free_imatrix(iyy, nrl, nrh, ncl, nch);
	free_imatrix(ixy, nrl, nrh, ncl, nch);
}

rgb8** display_corner(byte** I, int** C, long nrl, long nrh, long ncl, long nch)
{
	rgb8** disp_c = rgb8matrix(nrl, nrh, ncl, nch);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			
			disp_c[x][y].r = I[x][y]; disp_c[x][y].g = I[x][y]; disp_c[x][y].b = I[x][y];			
			if(C[x][y] == 1) {
				disp_c[x][y].r = 255;
				disp_c[x][y].g = 0; 
				disp_c[x][y].b = 0;				
			}
		}
	}		

	return disp_c;
}