// Libraries

#include <stdlib.h>
#include <stdio.h>

// Function declarations

int generate_problem(int, double, int *, int **, int **, double **, double **);
void generateOptions(int *, int *, int *, int *, int, char **);
void printAArrays(double *, int *, int *, int);
void printLinearSystemArrays(double *, int *, int *, double *, int);
int umfSolve(int, double *, int *j, int *, double *, double *);
