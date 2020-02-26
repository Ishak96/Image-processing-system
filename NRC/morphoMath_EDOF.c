#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

int ceil3(int x)
{
   return x * 255;
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

byte** erosion(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
	int same = 0;
	for (int i = 0; i < maskw; i++) {
		for (int j = 0; j < maskh; j++) {
			if(mask[i][j] == 1)
				same++;
		}
	}

   int** out = imatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         double acc = 0.0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
        
         out[x][y] = (acc == same) ? 1 : 0;
      }
   }

   apply(out, nrl, nrh, ncl, nch, ceil3);

   return convert_imatrix_bmatrix(out, nrl, nrh, ncl, nch);
}

byte** dilatation(byte** f, long nrl, long nrh, long ncl, long nch,
            float** mask, long maskw, long maskh)
{
   int** out = imatrix(nrl, nrh, ncl, nch);

   for (int x = nrl; x < nrh; x++) {
      for (int y = ncl; y < nch; y++) {
         double acc = 0.0;

         for (int u = 0; u < maskw; u++) {
            for (int v = 0; v < maskh; v++) {
               int nx = x + u - ((int) (maskw / 2));
               int ny = y + v - ((int) (maskh / 2));

               if(nx >= nrl && nx < nrh && ny >= ncl && ny < nch) {
                  acc += ((int)f[nx][ny] / 255) * mask[u][v];
               }

            }
         }
        
         out[x][y] = (acc > 0) ? 1 : 0;
      }
   }

   apply(out, nrl, nrh, ncl, nch, ceil3);

   return convert_imatrix_bmatrix(out, nrl, nrh, ncl, nch);
}

byte** minus(byte** I1, byte** I2,long nrl, long nrh, long ncl, long nch)
{
	int** out = imatrix(nrl, nrh, ncl, nch);

	for (int x = nrl; x < nrh; x++) {
		for (int y = ncl; y < nch; y++) {
			out[x][y] = (int)I1[x][y] - (int)I2[x][y];
		}
	}	

	return convert_imatrix_bmatrix(out, nrl, nrh, ncl, nch);
}

int main(void) {
	long nrh, nrl, nch, ncl;

	float maskE[3][3] = {{1, 1, 1},
      					  {1, 1, 1},
      					  {1, 1, 1}
   						 };
   	float** mask = imatrix(0, 3, 0, 3);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mask[i][j] = maskE[i][j];
		}
	}
	
	byte **I;
	int n = 6;

	byte** ITrou = LoadPGM_bmatrix("../Images/Test/carreTrou.pgm", &nrl, &nrh, &ncl, &nch);
	byte** IBruit = LoadPGM_bmatrix("../Images/Test/carreBruit.pgm", &nrl, &nrh, &ncl, &nch);

	//remplissage
	byte **IR = dilatation(ITrou, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		IR = dilatation(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		IR = erosion(IR, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	//debruitage
	n = 8;
	byte **ID = erosion(IBruit, nrl, nrh, ncl, nch, mask, 3, 3);
	for(int i = 1; i < n; i++) {
		ID = erosion(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}

	for(int i = 0; i < n; i++) {
		ID = dilatation(ID, nrl, nrh, ncl, nch, mask, 3, 3);
	}	

	SavePGM_bmatrix(IR, nrl, nrh, ncl, nch, "../Images/results/carre_Trou_remplis.pgm");
	SavePGM_bmatrix(ID, nrl, nrh, ncl, nch, "../Images/results/carre_Bruit_enleve.pgm");
	
	/* Free the byte sum matrix */
	free_bmatrix(ITrou, nrl, nrh, ncl, nch);
	free_bmatrix(IBruit, nrl, nrh, ncl, nch);
	free_bmatrix(IR, nrl, nrh, ncl, nch);
	free_bmatrix(ID, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, 3, 0, 3);

	return 0;
}