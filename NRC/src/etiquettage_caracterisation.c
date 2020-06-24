#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "etiquettage_caracterisation.h"
#include "utils.h"
#include "erosion_dilatation.h"
#include "extraction.h"
#include "labelling.h"

Parameters init_parameters()
{
	Parameters param = {
		.area_label = -1,
		.size = 0,
		.OX = 0, .OY = 0,
		.E_x = 0, .E_y = 0, .E_xy = 0,
		.theta = 0,
		.average_greyscal = 0, 
		.avrage_R = 0, .avrage_G = 0, .avrage_B = 0,
		.H = ivector(0, 255)
	};
	for (int i = 0; i < 255; i++)
		param.H[i] = 0;

	return param;
}

void free_parameters(Parameters *param)
{
	free_ivector(param->H, 0, 255);
}

void dump_parameters(Parameters* param, int size)
{
	printf("########IMAGE#########\n");
	for(int label = 0; label < size; label++) {
		printf("---area label--- : %d\n", param[label].area_label);
		printf("area size : %d\n", param[label].size);
		printf("mean x : %.2f\n", param[label].OX);
		printf("mean y : %.2f\n", param[label].OY);
		printf("mean E_x^2 : %.2f\n", param[label].E_x);
		printf("mean E_y^2 : %.2f\n", param[label].E_y);
		printf("mean E_xy^2 : %.2f\n", param[label].E_xy);
		printf("theta : %.2f\n", param[label].theta);
		printf("average greyscal : %.2f\n", param[label].average_greyscal);
		printf("average RGB [%.2f,%.2f,%.2f] :\n", param[label].avrage_R,
		 										   param[label].avrage_G, param[label].avrage_B);

		printf("[");
		for(int i = 0; i < 254; i++)
			printf("%d, ", param[label].H[i]);

		printf("%d]\n", param[label].H[254]);
	}
	printf("########IMAGE#########\n");
}

Parameters* compute_parameters(rgb8** IRGB, int** E, int num_label,
							  long nrl, long nrh, long ncl, long nch)
{
	Parameters* param = (Parameters *) malloc(num_label * sizeof(Parameters));

	for(int label = 1; label <= num_label; label++) {
		Parameters p = init_parameters();
		p.area_label = label;

		ArrayList x_list = NEW_ArrayList();
		ArrayList y_list = NEW_ArrayList();

		for (int y = nrl; y < nrh; y++) {
			for (int x = ncl; x < nch; x++) {
				if(E[y][x] == label) {
					p.size += 1;
					p.OX += x;
					p.OY += y;
					
					int byte = ((int)IRGB[y][x].r + (int)IRGB[y][x].g + (int)IRGB[y][x].b) / 3;
					p.average_greyscal += byte;

					p.avrage_R += (int)IRGB[y][x].r;
					p.avrage_G += (int)IRGB[y][x].g; 
					p.avrage_B += (int)IRGB[y][x].b;
					p.H[byte] += 1;

					ADD_ArrayList(&y_list, y); ADD_ArrayList(&x_list, x);
				}
			}
		}
		p.OX = p.OX / p.size;
		p.OY = p.OY / p.size;
		p.average_greyscal = p.average_greyscal / p.size;
		p.avrage_R = p.avrage_R / p.size;
		p.avrage_G = p.avrage_G / p.size;
		p.avrage_B = p.avrage_B / p.size;

		for(int i = 0; i < p.size; i++) {
			p.E_x += ((x_list.list[i] - p.OX) * (x_list.list[i] - p.OX)) / p.size;
			p.E_y += ((y_list.list[i] - p.OY) * (y_list.list[i] - p.OY)) / p.size;
			p.E_xy += ((x_list.list[i] - p.OX) * (y_list.list[i] - p.OY)) / p.size;
		}
		p.theta = 0.5 * atan((2 * p.E_x * p.E_y) / (p.E_x - p.E_y));

		param[label-1] = p;
		FREE_ArrayList(&x_list);
		FREE_ArrayList(&y_list);
	}

	return param;
}

void store_parameters(int i, int num_label, Parameters* param)
{
	const char* dir = "../Images/results/caracterisation/Lbox/";
	const char* label = "lbox";
	const char* extention = "txt";

	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);

	const char file_name[n];

	sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);

	FILE *fp;
	fp = fopen(file_name, "w");

	if(fp != NULL) {
		fprintf(fp, "NUMBER OF AREAS : %d\n", num_label);
		for(int label = 0; label < num_label; label++) {
			fprintf(fp, "---area label--- : %d\n", 			param[label].area_label);
			fprintf(fp, "area size : %d\n", 				param[label].size);
			fprintf(fp, "mean x : %.2f\n", 					param[label].OX);
			fprintf(fp, "mean y : %.2f\n", 					param[label].OY);
			fprintf(fp, "mean E_x^2 : %.2f\n", 				param[label].E_x);
			fprintf(fp, "mean E_y^2 : %.2f\n", 				param[label].E_y);
			fprintf(fp, "mean E_xy^2 : %.2f\n", 			param[label].E_xy);
			fprintf(fp, "theta : %.2f\n", 					param[label].theta);
			fprintf(fp, "average greyscal : %.2f\n", 		param[label].average_greyscal);
			fprintf(fp, "average RGB [%.2f,%.2f,%.2f]\n", param[label].avrage_R,
			      										    param[label].avrage_G, param[label].avrage_B);
			for(int i = 0; i < 254; i++) {
				if(param[label].H[i] != 0) {
					fprintf(fp,"(%d %d)", i, param[label].H[i]);
				}
			}

			fprintf(fp,"\n");
		}
	}

	fclose(fp);
}