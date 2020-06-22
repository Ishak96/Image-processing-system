#ifndef __CORNER_DETECTION_H__
#define __CORNER_DETECTION_H__

#define LAMBDA 0.005

rgb8** gradient(byte** I, long nrl, long nrh, long ncl, long nch);

rgb8** harris(byte** I, long nrl, long nrh, long ncl, long nch);

#endif