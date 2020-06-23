#ifndef __UTILS_H__
#define __UTILS_H__

#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

int pow2(int x);

int sqrt2(int x);

int I_min(int a, int b);

int I_max(int a, int b);

int ceil2(int x);

int ceil3(int x);

int** conv2(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh);

void apply(int** m, long nrl, long nrh, long ncl, long nch, int (*func)(int));

void bapply(byte** m, long nrl, long nrh, long ncl, long nch, int (*func)(int));

int** sum(int** I1, int** I2, long nrl, long nrh, long ncl, long nch);

void minus(byte** I1, byte** I2,long nrl, long nrh, long ncl, long nch);

int** multiply(int** I1, int** I2, long nrl, long nrh, long ncl, long nch);

void rgb8matrix_to_bmatrix(rgb8** I, byte** out, int nrl, int nrh, int ncl, int nch);

byte** convert_imatrix_bmatrix(int** m, long nrl, long nrh, long ncl, long nch);

byte** binarization(byte** I, long nrl, long nrh, long ncl, long nch);

int* histogram(byte** I, long nrl, long nrh, long ncl, long nch, long h_nrl, long h_nrh);

void rollingAverage(int* H, int h_nrl, int h_nrh, int N);

int** gaussian_filter(int** f, long nrl, long nrh, long ncl, long nch,
               		  float** mask, long maskw, long maskh);

int*** i3D(int M, int N, int O);

void free_i3D(int*** A, int M, int N);

int fileCount(const char* dir);

void print_progress(size_t count, size_t max);

#endif