#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "extraction.h"

/* ---------------------------------------- */
void erosion(byte** f, byte** out, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh)
/* ---------------------------------------- */
{
	int same = 0;
	for (int i = 0; i < maskw; i++) {
		for (int j = 0; j < maskh; j++) {
			if(mask[i][j] == 1)
				same++;
		}
	}

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         int acc = 0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
         out[x][y] = (acc == same) ? (byte) 255 : (byte) 0;
      }
   }
}

/* ---------------------------------------- */
void dilatation(byte** f, byte** out, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh)
/* ---------------------------------------- */
{
   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         int acc = 0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
         out[x][y] = (acc > 0) ? (byte) 255 : (byte) 0;
      }
   }
}

/* ---------------------------------------- */
void n_erosion(byte** src, byte** dist, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh, int n)
/* ---------------------------------------- */
{
   erosion(src, dist, nrl, nrh, ncl, nch, mask, maskw, maskh);
 
   for(int i = 1; i < n; i++) {
      copy_bmatrix(dist, nrl, nrh, ncl, nch, src, nrl, nrh, ncl, nch);
      erosion(src, dist, nrl, nrh, ncl, nch, mask, maskw, maskh);
   }
}

/* ---------------------------------------- */
void n_dilatation(byte** src, byte** dist, long nrl, long nrh, long ncl, long nch,
            int** mask, long maskw, long maskh, int n)
/* ---------------------------------------- */
{
   dilatation(src, dist, nrl, nrh, ncl, nch, mask, maskw, maskh);
 
   for(int i = 1; i < n; i++) {
      copy_bmatrix(dist, nrl, nrh, ncl, nch, src, nrl, nrh, ncl, nch);
      dilatation(src, dist, nrl, nrh, ncl, nch, mask, maskw, maskh);
   }
}