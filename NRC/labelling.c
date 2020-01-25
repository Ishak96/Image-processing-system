#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

#define INFINITY  10000

int CEIL = 110; // 110 || 35

/*mask functions to get label and attribut of A and B*/
int E_B(int** E, long nrl, long nrh, long ncl, long nch, int y, int x){
   int label = -1;

   if (y-1 > nrl) {
      label = E[y-1][x];
   }

   return label;
}

int E_A(int** E, long nrl, long nrh, long ncl, long nch, int y, int x){
   int label = -1;

   if (x-1 > ncl) {
      label = E[y][x-1];
   }

   return label;
}

int ATT_B(byte** I, long nrl, long nrh, long ncl, long nch, int y, int x){
   int att = -1;

   if (y-1 > nrl) {
      att = (int) I[y-1][x];
   }

   return att;
}

int ATT_A(byte** I, long nrl, long nrh, long ncl, long nch, int y, int x){
   int att = -1;

   if (x-1 > ncl) {
      att = (int) I[y][x-1];
   }

   return att;
}

/*intuitive algorithm of labelling*/
void CHANGE(int** E, long nrl, long nrh, long ncl, long nch, int e_a, int e_b)
{
   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         if (E[y][x] == e_a) {
            E[y][x] = e_b;
         }
      }
   }
}

int** INTUITIVE_LABELLING(byte** I, long nrl, long nrh, long ncl, long nch)
{
   int att_a, att_b, e_a, e_b, att_c;
   /*E : the matrix of labels*/
   int** E = imatrix(nrl, nrh, ncl, nch);

   for (int y = nrl; y < nrh; y++)
      for (int x = ncl; x < nch; x++)
         E[y][x] = 0;

   int CURRENTLABEL = 0;

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         /*recuperate the attribut and the label of A and B*/
         att_a = ATT_A(I, nrl, nrh, ncl, nch, y, x);
         att_b = ATT_B(I, nrl, nrh, ncl, nch, y, x);
         
         e_a = E_A(E, nrl, nrh, ncl, nch, y, x);
         e_b = E_B(E, nrl, nrh, ncl, nch, y, x);

         /*only if the attribut of C is different then 0 we'll apply the algorithm*/
         att_c = (int) I[y][x];
         if(att_c != 0){
            if (att_c == att_a && att_c != att_b){
               E[y][x] = e_a;
            } else if (att_c == att_b && att_c != att_a){
               E[y][x] = e_b;
            } else if (att_c != att_b && att_c != att_a) {
               CURRENTLABEL++;
               E[y][x] = CURRENTLABEL;
            } else if (att_c == att_b && att_c == att_a && e_a == e_b){
               E[y][x] = e_b;
            } else if (att_c == att_b && att_c == att_a && e_a != e_b) {
               E[y][x] = e_b;
               /*moves up in the picture and changes all the e_a to e_b*/
               CHANGE(E, nrl, nrh, ncl, nch, e_a, e_b);
            }
         }
      }
   }

   return E;
}

rgb8** CONVERT(int** E, long nrl, long nrh, long ncl, long nch)
{
   int R, G, B;
   int C;
   
   int a = (rand() % 256), b = (rand() % 256), c = (rand() % 256);

   rgb8** I = rgb8matrix(nrl, nrh, ncl, nch);

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {

         C = E[y][x];
         
         B = (C * a) % 256;
         G = (C * b) % 256;
         R = (C * c) % 256;
         
         I[y][x].r = R;
         I[y][x].g = G;
         I[y][x].b = B;
      }
   }

   return I;
}

rgb8** EXTRACTZONE(byte** I, int** E, int label, long nrl, long nrh, long ncl, long nch)
{
   rgb8** R = rgb8matrix(nrl, nrh, ncl, nch);

   int dx = 0, dy = 0, N = 0;
   int min_x = INFINITY, min_y = INFINITY, max_x = -INFINITY, max_y = -INFINITY;

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
            R[y][x].r = I[y][x];
            R[y][x].g = I[y][x];
            R[y][x].b = I[y][x];
      }
   }

   for (int y = nrl; y < nrh; y++) {
      for (int x = ncl; x < nch; x++) {
         if (E[y][x] == label) {
            dx += x;
            dy += y;
            N++;
            
            if(x > max_x) max_x = x;
            if(x < min_x) min_x = x;
            
            if(y > max_y) max_y = y;
            if(y < min_y) min_y = y;
         }
      }
   }

   dx /= N; dy /= N;

   R[dy][dx].r = 255; R[dy][dx].b = 0; R[dy][dx].b = 0;

   for (int y = min_y; y < max_y; y++) {
      for (int x = min_x; x < max_x; x++) {
         if(x >= max_x - 1 || y >= max_y - 1 || x <= min_x  || y <= min_y) {
            R[y][x].r = 255; R[y][x].b = 0; R[y][x].b = 0;
         }
      }
   }

   return R;
}

/*algorithm with look-up table*/

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

   I = LoadPGM_bmatrix("../Images/src/rice.pgm", &nrl, &nrh, &ncl, &nch);

   byte** R = binarization(I, nrl, nrh, ncl, nch);

   int** E = INTUITIVE_LABELLING(R, nrl, nrh, ncl, nch);

   rgb8** L = CONVERT(E, nrl, nrh, ncl, nch);
   
   SavePPM_rgb8matrix(L, nrl, nrh, ncl, nch, "../Images/results/rice_labelling_intuitive.ppm");

   int label = 299;

   rgb8** EX = EXTRACTZONE(I, E, label, nrl, nrh, ncl, nch);
   SavePPM_rgb8matrix(EX, nrl, nrh, ncl, nch, "../Images/results/rice_zone20.pgm");


   /* Free the byte sum matrix */
   free_bmatrix(I, nrl, nrh, ncl, nch);
   free_bmatrix(R, nrl, nrh, ncl, nch);
   free_imatrix(E, nrl, nrh, ncl, nch);
   
   free_rgb8matrix(L, nrl, nrh, ncl, nch);
   free_rgb8matrix(EX, nrl, nrh, ncl, nch);

   return 0;
}