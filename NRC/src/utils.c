#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"

int CEIL;

int pow2(int x)
{
	return x * x;
}

int sqrt2(int x)
{
	return (int) sqrt((double) x);
}

int I_min(int a, int b)
{
	return (a > b) ? b : a;
}

int I_max(int a, int b)
{
	return (a > b) ? a : b;
}

int ceil2(int x)
{
	return 255 * (x > CEIL) ;
}

int ceil3(int x)
{
   return x * 255;
}

/* 2D convolution with a mask
 * Retuns an imatrix to avoid overflows
 */
int** conv2(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
	int** out = imatrix(nrl, nrh, ncl, nch);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			double acc = 0.0;
			int n = 0;

			for (int u = 0; u < maskw; u++) {
				for (int v = 0; v < maskh; v++) {
					int nx = x + u - ((int) (maskw / 2));
					int ny = y + v - ((int) (maskh / 2));

					if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
						acc += f[nx][ny] * mask[u][v];
						n++;
					}
				}
			}
			out[x][y] = acc / n;
		}
	}
	
	return out;
}

/* Apply a function to every element of an imatrix*/
void apply(int** m, long nrl, long nrh, long ncl, long nch, int (*func)(int))
{
	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			m[x][y] = func(m[x][y]);
		}
	}
}

/* Apply a function to every element of an byte matrix*/
void bapply(byte** m, long nrl, long nrh, long ncl, long nch, int (*func)(int))
{

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         m[y][x] = func(m[y][x]);
      }
   }
}

int** sum(int** I1, int** I2, long nrl, long nrh, long ncl, long nch)
{
	int** result = imatrix(nrl, nrh, ncl, nch);
	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			result[x][y] = I1[x][y] + I2[x][y];
		}
	}
	return result;
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

int** multiply(int** I1, int** I2, long nrl, long nrh, long ncl, long nch)
{
	int** result = imatrix(nrl, nrh, ncl, nch);
	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			int res = I1[x][y] * I2[x][y];
			result[x][y] = res;
		}
	}
	return result;
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

/* Convert an imatrix to a bmatrix
 * PS: casting to smaller type can be lossy
 */
byte** convert_imatrix_bmatrix(int** m, long nrl, long nrh, long ncl, long nch)
{

   byte** out = bmatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         out[x][y] = (byte) m[x][y];
      }
   }

   free_imatrix(m, nrl, nrh, ncl, nch);

   return out;
}

byte** binarization(byte** I, long nrl, long nrh, long ncl, long nch)
{
   byte** out = bmatrix(nrl, nrh, ncl, nch);
   
   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         out[y][x] = I[y][x];
      }
   }  

   bapply(out, nrl, nrh, ncl, nch, ceil2);

   return out;
}

int* histogram(byte** I, long nrl, long nrh, long ncl, long nch, long h_nrl, long h_nrh)
{
	int* H = ivector(h_nrl, h_nrh);

	int index;

	//initialization of the vector
	for (int x = h_nrl; x < h_nrh; x++) {
		H[x] = 0;
	}

	//histogram process
	for (int x = nrl + 1; x < nrh - 1; x++) {
		for (int y = ncl + 1; y < nch - 1; y++) {
			index = I[x][y];
			H[index]++;
		}
	}  
	return H;
}

void rollingAverage(int* H, int h_nrl, int h_nrh, int N)
{
	int som;

	if(N > 0)
		for (int x = h_nrh -1; x >= h_nrl; x--) {
			if (x >= (h_nrl + N)) {
				som = 0;
				for (int t = 0; t < N; t++) {
					som += H[x-t];
				}
				H[x] = som / N;
			}
		}
}

int** gaussian_filter(int** f, long nrl, long nrh, long ncl, long nch,
               		  float** mask, long maskw, long maskh)
{
	int** out = imatrix(nrl, nrh, ncl, nch);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			double acc = 0.0;
			float n = 0.0;

			for (int u = 0; u < maskw; u++) {
				for (int v = 0; v < maskh; v++) {
					int nx = x + u - ((int) (maskw / 2));
					int ny = y + v - ((int) (maskh / 2));

					if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
						acc += f[nx][ny] * mask[u][v];
						n += mask[u][v];
					}
				}
			}
			out[x][y] = (acc / 16);
		}
	}
	return out;
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

int fileCount(const char* dir) 
{
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