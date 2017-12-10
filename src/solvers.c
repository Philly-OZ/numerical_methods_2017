// Author : Denis Verstraeten
// Created on : 08/12/2017

/* This module contains all the solving methods */

#include "solvers.h"

int umfSolve(int problemSize, double *a, int *ja, int *ia, double *x,
  double *b){
    // this is called to solve the linear problem Ax=b
    void *Numeric = NULL; // variable required for the factorisation
    if (factorisation(problemSize, a, ja, ia, &Numeric) ||
  resolution(ia, ja, a, x, b, &Numeric)){
      return EXIT_FAILURE;
    } else {
      return EXIT_SUCCESS;
    }
  }

int sgsSolve(double *a, int *ia, int *ja, double **x, double *b,
  double minResidue, int maxIterations, int m, int problemSize, int PREC_DEBUG){
    /* this function solves Ax=b iteratively using Symmetric Gauss Seidel
    preconditionning.
    a, ia, ja : arrays containing the A matrix in CSR format
    x : vector of solution
    b : vector of independent terms
    minResidue : tolerance on the residue, the algorithm stops when the residue
    is smaller than minResidue
    maxIterations : maximum number of iteraions allowed (to avoid infinite
    loops if the algorithm diverges)
    m : discretisation of the problem
    problemSize : number of unknowns
    PREC_DEBUG : binary variable indicating if code is in debug about the
    preconditionning */

    // Generating preconditionner B^-1 = U^-1 D L^-1

    // First, generating U, D and L

    double *la, *ua, *da; /* arrays that will be filled by the U, L and D in CSR
    format */
    int *ila, *jla, *iua, *jua; /* arrays that will be filled by the U, L and D
    in CSR format */

    if(splitAMatrix(m, problemSize, a, ia, ja, &la, &ila, &jla, &ua, &iua, &jua,
    &da)){
      printf("ERROR : Splitting matrix failed.\n");
      return EXIT_FAILURE;
    }

    if (PREC_DEBUG){
      // code is in preconditionner debugging mode for SGS solving
      printf("DEBUGGING : preconditionning for SGS\n");
      printf("Printing the arrays...\n\n");
      printPrecSGSArrays(la, ua, ila, iua, jla, jua, da, problemSize);
    }

    return EXIT_SUCCESS;
  }
