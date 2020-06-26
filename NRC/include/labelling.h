#ifndef __LABELLING_H__
#define __LABELLING_H__

typedef int TElement;

#define INFINITY  10000

typedef struct
{
   int size;
   int current_size;
   TElement* list;
}ArrayList;

/*functions for ArrayList*/
ArrayList NEW_ArrayList();

void ADD_ArrayList(ArrayList* arr, TElement elem);

void FREE_ArrayList(ArrayList* arr);

/*mask functions to get label and attribut of A and B*/
int E_B(int** E, long nrl, long nrh, long ncl, long nch, int y, int x);

int E_A(int** E, long nrl, long nrh, long ncl, long nch, int y, int x);

int ATT_B(byte** I, long nrl, long nrh, long ncl, long nch, int y, int x);

int ATT_A(byte** I, long nrl, long nrh, long ncl, long nch, int y, int x);

/*intuitive algorithm of labelling*/
void CHANGE(int** E, long nrl, long nrh, long ncl, long nch, int e_a, int e_b);

int** INTUITIVE_LABELLING(byte** I, int* CURRENTLABEL, long nrl, long nrh, long ncl, long nch);

rgb8** CONVERT(int** E, long nrl, long nrh, long ncl, long nch);

void EXTRACTZONE(rgb8** R, byte** I, int** E, int label, long nrl, long nrh, long ncl, long nch);

void EXTRACTZONE_RGB8(rgb8** R, int** E, int label, long nrl, long nrh, long ncl, long nch);

/*algorithm with look-up table*/
void LOOKUP_TABLE_LABELLING(byte** I, int** E, int* CURRENTLABEL, long nrl, long nrh, long ncl, long nch);

#endif