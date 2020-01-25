#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

#define INFINITY  10000

int CEIL = 110;

/* Apply a function to every element of an imatrix
 */
void apply(byte** m, long nrl, long nrh, long ncl, long nch, int (*func)(int))
{

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         m[x][y] = func(m[x][y]);
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
	
   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         out[x][y] = I[x][y];
      }
   }  

   apply(out, nrl, nrh, ncl, nch, ceil2);

   return out;
}

int* histogram(byte** I, long nrl, long nrh, long ncl, long nch, long h_nrl, long h_nrh)
{
   int* H = ivector(h_nrl, h_nrh);

   int index;

   //initialization of the vector
   for (int x = h_nrl; x < h_nrh; x++) {
      H[x] = 0;
   }

   //histogram process
   for (int x = nrl + 1; x < nrh - 1; x++) {
      for (int y = ncl + 1; y < nch - 1; y++) {
         index = I[x][y];
         H[index]++;
      }
   }  

   return H;
}

void rollingAverage(int* H, int h_nrl, int h_nrh, int N)
{
   int som;

   if(N > 0)
      for (int x = h_nrh -1; x >= h_nrl; x--) {
         if (x >= (h_nrl + N)) {
            som = 0;
            for (int t = 0; t < N; t++) {
               som += H[x-t];
            }
            H[x] = som / N;
         }
      }
}

int computeCeil(int* H, int h_nrl, int h_nrh)
{
   int max1 = -INFINITY, max2 = -INFINITY;
   int index1 = 0, index2 = 0;

   for (int x = h_nrl; x < h_nrh; x++) {

      if (H[x] > max1) {
         max1 = H[x];
         index1 = x;

         index2 = index1;
      } else if (H[x] > max2 && x > index1+10 && H[x] > 500) {
         max2 = H[x];
         index2 = x;
      }
   }

   return (H[index1] > H[index2]) ? index2 : index1;;
}

int main(void) {
   long nrh, nrl,
        nch, ncl;

   long h_nrl = 0;
   long h_nrh = 255;

   int N = 15;

   byte **I;

   //cubesx3 || rice
   I = LoadPGM_bmatrix("../Images/src/rice.pgm", &nrl, &nrh, &ncl, &nch);

   int* H = histogram(I, nrl, nrh, ncl, nch, h_nrl, h_nrh);

   rollingAverage(H, h_nrl, h_nrh, N);

   CEIL = computeCeil(H, h_nrl, h_nrh);

   byte** R = binarization(I, nrl, nrh, ncl, nch);

   SavePGM_bmatrix(R, nrl, nrh, ncl, nch, "../Images/results/rice_binarized.pgm");

   /* Free the byte sum matrix */
   free_bmatrix(I, nrl, nrh, ncl, nch);
   free_bmatrix(R, nrl, nrh, ncl, nch);
   free_ivector(H, h_nrl, h_nrh);

   return 0;
}