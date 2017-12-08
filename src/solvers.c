// Author : Denis Verstraeten
// Created on : 08/12/2017

/* This module contains all the solving methods */

#include "solvers.h"

int umfSolve(int problemSize, double *a, int *ja, int *ia, double *x, \
  double *b){
    // this is called to solve the linear problem Ax=b
    void *Numeric = NULL; // variable required for the factorisation
    if (factorisation(problemSize, a, ja, ia, &Numeric) || \
  resolution(ia, ja, a, x, b, &Numeric)){
      return EXIT_FAILURE;
    } else {
      return EXIT_SUCCESS;
    }
  }
