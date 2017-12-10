// Libraries

#include <stdlib.h>
#include <stdio.h>
#include "umfpack.h"

// Defines

#define square(x) ((x) * (x))

// Function declaration

double sqrt(double);
double exp(double);
int umfSolve(int, double *, int *, int *, double *, double *);
