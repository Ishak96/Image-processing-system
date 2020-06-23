#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "labelling.h"
#include "utils.h"

extern int CEIL = 110; // 110 || 35

int main(void) {
   long nrh, nrl,
        nch, ncl;

   byte **I;
   int nb_label;

   I = LoadPGM_bmatrix("../Images/src/rice.pgm", &nrl, &nrh, &ncl, &nch);

   int** E = imatrix(nrl, nrh, ncl, nch);

   byte** R = binarization(I, nrl, nrh, ncl, nch);

   LOOKUP_TABLE_LABELLING(R, E, &nb_label, nrl, nrh, ncl, nch);

   rgb8** L = CONVERT(E, nrl, nrh, ncl, nch);
   
   SavePPM_rgb8matrix(L, nrl, nrh, ncl, nch, "../Images/results/rice_labelling_lookup.ppm");

   int label = 7;

   rgb8** EX = EXTRACTZONE(I, E, label, nrl, nrh, ncl, nch);
   SavePPM_rgb8matrix(EX, nrl, nrh, ncl, nch, "../Images/results/rice_zone7.pgm");

   /* Free the byte sum matrix */
   free_bmatrix(I, nrl, nrh, ncl, nch);
   free_bmatrix(R, nrl, nrh, ncl, nch);
   free_imatrix(E, nrl, nrh, ncl, nch);
   
   free_rgb8matrix(L, nrl, nrh, ncl, nch);
   free_rgb8matrix(EX, nrl, nrh, ncl, nch);

   return 0;
}