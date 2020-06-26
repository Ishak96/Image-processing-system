#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "extraction.h"

/* ---------------------------------------- */
byte** temporalAverage(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq)
/* ---------------------------------------- */
{
	int n = strlen(dir_name) + strlen(label) + 3 + strlen(extention);
	
	double** out = dmatrix(*nrl, *nrh, *ncl, *nch);;

	const char file_name[n];

	for(int i = 1; i <= n_seq ; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir_name, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, nrl, nrh, ncl, nch);

		for (int x = *nrl; x < *nrh; x++) {
			for (int y = *ncl; y < *nch; y++) {
				double value = ((int)I1[x][y].r + (int)I1[x][y].g + (int)I1[x][y].b) / 3;
				out[x][y] = out[x][y] + (value / n_seq);
			}
		}

		free_rgb8matrix(I1, *nrl, *nrh, *ncl, *nch);
	}

	byte** bOut = bmatrix(*nrl, *nrh, *ncl, *nch);			

	for (int x = *nrl; x < *nrh; x++) {
		for (int y = *ncl; y < *nch; y++) {
			bOut[x][y] = (byte) abs(out[x][y]);
		}
	}

	free_dmatrix(out, *nrl, *nrh, *ncl, *nch);

	return bOut;
}

/* ---------------------------------------- */
byte** medianFilter(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq)
/* ---------------------------------------- */
{
	int n = strlen(dir_name) + strlen(label) + 3 + strlen(extention);
	const char file_name[n];

	int som, median;
	
	int*** votingBox;

	for(int i = 1; i <= n_seq ; i++) {
		sprintf(file_name, "%s%s%03d.%s", dir_name, label, i, extention);
		rgb8** I1 = LoadPPM_rgb8matrix(file_name, nrl, nrh, ncl, nch);

		if(i == 1)
			votingBox = i3D(*nrh, *nch, 256);			
		
		for (int x = *nrl; x < *nrh; x++) {
			for (int y = *ncl; y < *nch; y++) {
				int value = ((int)I1[x][y].r + (int)I1[x][y].g + (int)I1[x][y].b) / 3;
				votingBox[x][y][value]++;
			}
		}

		free_rgb8matrix(I1, *nrl, *nrh, *ncl, *nch);
	}

	byte** bOut = bmatrix(*nrl, *nrh, *ncl, *nch);

	median = n_seq / 2;
	int find = 0;

	for (int x = *nrl; x < *nrh; x++) {
		for (int y = *ncl; y < *nch; y++) {
			som = 0;
			find = 0;
			for(int z = 0; z < 256 && !find; z++) {
				som += votingBox[x][y][z];

				if(som > median) {
					bOut[x][y] = (byte) z;
					find = 1;
				}
			}
		}
	}

	free_i3D(votingBox, *nrh, *nch);

	return bOut;
}