#ifndef __EROSION_DILATATION_H__
#define __EROSION_DILATATION_H__


byte** erosion(byte** f, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh);

byte** dilatation(byte** f, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh);

#endif