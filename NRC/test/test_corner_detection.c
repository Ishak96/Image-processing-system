#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "corner_detection.h"
#include "utils.h"

int main(void) {
	long nrh, nrl,
         nch, ncl;
	byte **I;

	I = LoadPGM_bmatrix("../Images/src/cubesx3.pgm", &nrl, &nrh, &ncl, &nch);

	rgb8** harris_i = harris(I, nrl, nrh, ncl, nch);

	SavePPM_rgb8matrix(harris_i, nrl, nrh, ncl, nch, "../Images/results/corner_detection/gradient_corner_detection_cubesx3.ppm");

	/* Free the byte sum matrix */
	free_bmatrix(I, nrl, nrh, ncl, nch);
	free_rgb8matrix(harris_i, nrl, nrh, ncl, nch);
	
	return 0;
}