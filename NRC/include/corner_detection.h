#ifndef __CORNER_DETECTION_H__
#define __CORNER_DETECTION_H__

#define LAMBDA 0.00005

int** gradient(byte** I, long nrl, long nrh, long ncl, long nch);

int** harris(byte** I, long nrl, long nrh, long ncl, long nch);

rgb8** display_corner(byte** I, int** C, long nrl, long nrh, long ncl, long nch);

#endif