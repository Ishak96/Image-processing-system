#ifndef __CORNER_DETECTION_H__
#define __CORNER_DETECTION_H__

#define LAMBDA 0.00005

void create_mask(int harris_gradient);

void destroy_mask();

int** gradient(byte** I, long nrl, long nrh, long ncl, long nch);

void harris(int** harris_i, byte** I, long nrl, long nrh, long ncl, long nch);

rgb8** display_corner(byte** I, int** C, long nrl, long nrh, long ncl, long nch);

#endif