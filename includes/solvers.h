// Libraries

#include <stdlib.h>
#include <stdio.h>

// Function declaration

int factorisation(int, double *, int *, int *, void **);
int resolution(int *, int *, double *, double *, double *, void **);
int splitAMatrix(int, int, double *, int *, int *, double **, int **, int **,
  double **, int **, int **, double **, int **, int **);
void printSplitSGSArrays(double *, double *, int *, int *, int *, int *,
  double *, int *, int *, int);
int inverseMatrix(int, double **, int **, int **, double *, int *, int *, int);
void printInvSGSArrays(double *, double *, int *, int *, int *, int *, int);
void printPrecSGSArrays(double *, int *, int *, int);
int getResidue(int, double *, double *, int *, int *, double *, double **);
double normVector(int , double *);
int matrixVectorMultCSR(int, double *, int *, int *, double *, double **);
int vectorAddition(int, double *, double *, double **);
int saveResidue(FILE *, int, double);
int plotResidue(char *, char *);
