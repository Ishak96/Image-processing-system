#ifndef __TRACK_H__
#define __TRACK_H__

#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "etiquettage_caracterisation.h"

typedef struct {
	
	int OX, OY;
	int min_x, max_x; 
	int min_y, max_y;
} Pose;

typedef struct {
	
	int size;
	Pose pose;
	int area_size;
	int* H;
} Vignette;

void init_system(Vignette* vignette, rgb8** I, int** harris_i, int** E, 
				 int track_label, int vignette_size, long nrl, long nrh, long ncl, long nch);

void compute_param(rgb8** I, int** E, int label, int* H, int* size,
				   int* min_x, int* max_x, int* min_y, int* max_y,
				   long nrl, long nrh, long ncl, long nch);

void dump_vignette(Vignette vignette);

int track(int** E, rgb8** I1, Vignette* vignette, long nrl, long nrh, long ncl, long nch);

void extract_zone_p(rgb8** I, byte** byteI1, Parameters param,
				  int rgb, long nrl, long nrh, long ncl, long nch);

void extract_zone(rgb8** I, byte** byteI1, Vignette vignette,
				  int rgb, long nrl, long nrh, long ncl, long nch);

#endif