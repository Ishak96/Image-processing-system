#include "track.h"

/* ---------------------------------------- */
void init_system(Vignette* vignette, rgb8** I, int** harris_i, int** E, 
				 int track_label, int vignette_size, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------- */
{
	vignette->pose.OY = -INFINITY;
	int N = 0;

	vignette->H = ivector(0, 255);
	for (int i = 0; i < 255; i++)
		vignette->H[i] = 0;	

	for (int y = nrl; y < nrh; y++) {
		for (int x = ncl; x < nch; x++) {
			if (E[y][x] == track_label) {
				int value = ((int)I[x][y].r + (int)I[x][y].g + (int)I[x][y].b) / 3;
            	vignette->H[value] += 1;
            	N++;
            
				if(y > vignette->pose.OY && harris_i[y][x]) {
					vignette->pose.OX = x;
					vignette->pose.OY = y;
				}
			}
		}
	}

	vignette->size = vignette_size;
	vignette->area_size = N;
	vignette->pose.min_x = vignette->pose.max_x = 0;
	vignette->pose.min_y = vignette->pose.max_y = 0;
}

/* ---------------------------------------- */
void compute_param(rgb8** I, int** E, int label, int* H, int* size,
				   int* min_x, int* max_x, int* min_y, int* max_y,
				   long nrl, long nrh, long ncl, long nch) 
/* ---------------------------------------- */
{	

	for (int i = 0; i < 255; i++)
		H[i] = 0;

	*size = 0;

	for (int y = nrl; y < nrh; y++) {
		for (int x = ncl; x < nch; x++) {
			if (E[y][x] == label) {
				
				if(x > *max_x) *max_x = x;
				if(x < *min_x) *min_x = x;
            
				if(y > *max_y) *max_y = y;
				if(y < *min_y) *min_y = y;

				*size = *size + 1;
				
				int value = ((int)I[y][x].r + (int)I[y][x].g + (int)I[y][x].b) / 3;
            	H[value] += 1;

            }
        }
    }
}

/* ---------------------------------------- */
void dump_vignette(Vignette vignette)
/* ---------------------------------------- */
{
	printf("size : %d\n", vignette.size);
	for(int i = 0; i < 255; i++)
		printf("[%d->%d]", i, vignette.H[i]);
	printf("\n");
}

/* ---------------------------------------- */
int track(int** E, rgb8** I1, Vignette* vignette, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------- */
{
	int new_label = 0;
	int h_err = 10000;
	int app_size = 10000;

	int new_x = vignette->pose.OX, new_y = vignette->pose.OY;
	
	int curr_min_x = vignette->pose.min_x;
	int curr_max_x = vignette->pose.max_x;
	int curr_min_y = vignette->pose.min_y;
	int curr_max_y = vignette->pose.max_y;
	int curr_area_size = vignette->area_size;

	int* H = ivector(0, 255);

	for (int u = 0; u < vignette->size; u++) {
		for (int v = 0; v < vignette->size; v++) {
				
			// The thumbnail centred on the point of interest
			int ny = vignette->pose.OY + u - ((int) (vignette->size / 2));
			int nx = vignette->pose.OX + v - ((int) (vignette->size / 2));

			if(ny >= nrl && ny < nrh && nx >= ncl && nx < nch) {
					
				int label = E[ny][nx];
				if(label != 0) {
					int next_min_x = INFINITY, next_max_x = -INFINITY, 
						next_min_y = INFINITY, next_max_y = -INFINITY;
					int next_area_size = 0;

					compute_param(I1, E, label, H, &next_area_size,
								  &next_min_x, &next_max_x, &next_min_y, &next_max_y,
								  nrl, nrh, ncl, nch);
		
					int curr_h_err = 0;
					int curr_app_size = 100 * (abs(next_area_size - curr_area_size) 
											/ (float)curr_area_size);
					
					for(int j = 0; j < 255; j++) {
						if(vignette->H[j] != 0 && H[j] != vignette->H[j]) {
							curr_h_err++;
						}
					}

					if(curr_h_err < h_err && curr_app_size < app_size) {
						app_size = curr_app_size;
						h_err = curr_h_err;
						new_label = label;
						
						new_y = ny;
						new_x = nx;
						
						curr_min_x = next_min_x;
						curr_max_x = next_max_x;
						curr_min_y = next_min_y;
						curr_max_y = next_max_y;
					}
				}
			}
		}
	}
	
	free_ivector(H, 0, 255);

	vignette->pose.OY = new_y;
	vignette->pose.OX = new_x;

	vignette->pose.min_x = curr_min_x;
	vignette->pose.max_x = curr_max_x;
	vignette->pose.min_y = curr_min_y;
	vignette->pose.max_y = curr_max_y;	

	return new_label;
}

/* ---------------------------------------- */
void extract_zone_p(rgb8** I, byte** byteI1, Parameters param,
				  int rgb, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------- */
{
	int min_x = param.min_x, min_y = param.min_y; 
	int max_x = param.max_x, max_y = param.max_y;
	
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
				I[y][x].r = 0; I[y][x].g = 255; I[y][x].b = 0;
			}
		}
	}
}

/* ---------------------------------------- */
void extract_zone(rgb8** I, byte** byteI1, Vignette vignette,
				  int rgb, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------- */
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
	I[vignette.pose.OY][vignette.pose.OX].r = 0; 
	I[vignette.pose.OY][vignette.pose.OX].g = 0; 
	I[vignette.pose.OY][vignette.pose.OX].b = 255;
	for (int y = min_y; y < max_y; y++) {
		for (int x = min_x; x < max_x; x++) {
			if(x >= max_x - 1 || y >= max_y - 1 || x <= min_x  || y <= min_y) {
				I[y][x].r = 255; I[y][x].g = 0; I[y][x].b = 0;
			}
		}
	}
}