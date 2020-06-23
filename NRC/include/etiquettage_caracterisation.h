#ifndef __E_C_H__
#define __E_C_H__

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

Parameters init_parameters();

void free_parameters(Parameters *param);

void dump_parameters(Parameters* param, int size);

Parameters* compute_parameters(rgb8** IRGB, int** E, int num_label,
							  long nrl, long nrh, long ncl, long nch);

void store_parameters(int i, int num_label, Parameters* param);

#endif