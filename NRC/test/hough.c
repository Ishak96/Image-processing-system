#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

int CEIL = 110;

/* Apply a function to every element of an imatrix*/
void apply(byte** m, long nrl, long nrh, long ncl, long nch, int (*func)(int))
{

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         m[y][x] = func(m[y][x]);
      }
   }
}

int ceil2(int x)
{
   return 255 * (x > CEIL) ;
}

byte** binarization(byte** I, long nrl, long nrh, long ncl, long nch)
{
   byte** out = bmatrix(nrl, nrh, ncl, nch);
   
   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         out[y][x] = I[y][x];
      }
   }  

   apply(out, nrl, nrh, ncl, nch, ceil2);

   return out;
}

int main(void) {
   long nrh, nrl,
        nch, ncl;

   byte **I;

   I = LoadPGM_bmatrix("../Images/src/route.pgm", &nrl, &nrh, &ncl, &nch);

   byte** R = binarization(I, nrl, nrh, ncl, nch);

   /* Free the byte sum matrix */
   free_bmatrix(I, nrl, nrh, ncl, nch);
   free_bmatrix(R, nrl, nrh, ncl, nch);
   
   return 0;
}