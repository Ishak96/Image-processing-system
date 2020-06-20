#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"


/* 2D convolution with a mask
 * Retuns an imatrix to avoid overflows
 */
int** conv2(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
   int** out = imatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         double acc = 0.0;
         int n = 0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += f[nx][ny] * mask[u][v];
                  n++;
               }

            }
         }

         out[x][y] = acc / n;
      }
   }

   return out;
}

/* Apply a function to every element of an imatrix
 */
void apply(int** m, long nrl, long nrh, long ncl, long nch, int (*func)(int))
{

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         m[x][y] = func(m[x][y]);
      }
   }
}

/* Add two imatrices
 */
int** add(int **m1, int **m2, long nrl, long nrh, long ncl, long nch)
{
   int** out = imatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         out[x][y] = m1[x][y] + m2[x][y];
      }
   }

   return out;
}

/* Convert an imatrix to a bmatrix
 * PS: casting to smaller type can be lossy
 */
byte** convert_imatrix_bmatrix(int** m, long nrl, long nrh, long ncl, long nch)
{

   byte** out = bmatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         out[x][y] = (byte) m[x][y];
      }
   }

   free_imatrix(m, nrl, nrh, ncl, nch);

   return out;
}

int pow2(int x)
{
   return x * x;
}

int sqrt2(int x)
{
   return (int) sqrt((double) x);
}

int ceil2(int x)
{
   return 255 * (x > 5) ;
}

byte** bgradient(byte** I, long nrl, long nrh, long ncl, long nch)
{
   float maskd1[3][3] = {{-1, -2, -1},
      {0, 0, 0},
      {1, 2, 1}
   };
   float maskd2[3][3] = {{-1, 0, 1},
      {-2, 0, 2},
      {-1, 0, 1}
   };

   float** mask1 = malloc(sizeof(float*) * 3);
   float** mask2 = malloc(sizeof(float*) * 3);

   for (int i = 0; i < 3; i++) {
      mask1[i] = malloc(sizeof(float) * 3);
      mask2[i] = malloc(sizeof(float) * 3);

      for (int j = 0; j < 3; j++) {
         mask1[i][j] = maskd1[i][j];
         mask2[i][j] = maskd2[i][j];
      }
   }

   /* Calculate convolutions for the two masks */
   int** im1 = conv2(I, nrl, nrh, ncl, nch, mask1, 3, 3); // horizontal
   int** im2 = conv2(I, nrl, nrh, ncl, nch, mask2, 3, 3); // vertical

   /* Normalize */
   /* apply square */
   apply(im1, nrl, nrh, ncl, nch, pow2);
   apply(im2, nrl, nrh, ncl, nch, pow2);

   /* add two matrices */
   int** sum = add(im1, im2, nrl, nrh, ncl, nch);
   /* apply square root */
   apply(sum, nrl, nrh, ncl, nch, sqrt2);

   /* Apply ceil */
   apply(sum, nrl, nrh, ncl, nch, ceil2);

   /* Free the dynamically allocated masks */
   for (int i = 0; i < 3; i++) {
      free(mask1[i]);
      free(mask2[i]);
   }

   free(mask1);
   free(mask2);

   /* Free the convolution results*/
   free_imatrix(im1, nrl, nrh, ncl, nch);
   free_imatrix(im2, nrl, nrh, ncl, nch);

   /* convert to byte matrix and save (and frees the sum matrix) */
   return convert_imatrix_bmatrix(sum, nrl, nrh, ncl, nch);
}

rgb8** rgb8gradient(rgb8** m, long nrl, long nrh, long ncl, long nch)
{
   byte **R = bmatrix(nrl, nrh, ncl, nch),
        **G = bmatrix(nrl, nrh, ncl, nch),
        **B = bmatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         R[x][y] = m[x][y].r;
         G[x][y] = m[x][y].g;
         B[x][y] = m[x][y].b;
      }
   }

   byte** R_Grad = bgradient(R, nrl, nrh, ncl, nch);
   byte** G_Grad = bgradient(G, nrl, nrh, ncl, nch);
   byte** B_Grad = bgradient(B, nrl, nrh, ncl, nch);

   free_bmatrix(R, nrl, nrh, ncl, nch);
   free_bmatrix(G, nrl, nrh, ncl, nch);
   free_bmatrix(B, nrl, nrh, ncl, nch);


   rgb8** out = rgb8matrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         out[x][y].r = out[x][y].g = out[x][y].b = 255 * ((R_Grad[x][y] == 255) + (G_Grad[x][y] == 255) + (B_Grad[x][y] == 255) > 2);
      }
   }

   free_bmatrix(R_Grad, nrl, nrh, ncl, nch);
   free_bmatrix(G_Grad, nrl, nrh, ncl, nch);
   free_bmatrix(B_Grad, nrl, nrh, ncl, nch);

   return out;
}

int main(void) {
   long nrh, nrl,
        nch, ncl;
   rgb8 **I;

   I = LoadPPM_rgb8matrix("../Images/ppm/kids.ppm", &nrl, &nrh, &ncl, &nch);

   rgb8** Grad = rgb8gradient(I, nrl, nrh, ncl, nch);
   SavePPM_rgb8matrix(Grad, nrl, nrh, ncl, nch, "../Images/results/kids_gredient.ppm");

   /* Free the byte sum matrix */
   free_rgb8matrix(Grad, nrl, nrh, ncl, nch);
   free_rgb8matrix(I, nrl, nrh, ncl, nch);

   return 0;
}
