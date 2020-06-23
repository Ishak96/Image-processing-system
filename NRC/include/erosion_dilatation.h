#ifndef __EROSION_DILATATION_H__
#define __EROSION_DILATATION_H__


void erosion(byte** f, byte** out, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh);

void dilatation(byte** f, byte** out, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh);

void n_erosion(byte** src, byte** dist, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh, int n);

void n_dilatation(byte** src, byte** dist, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh, int n);

#endif