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

    double *la, *ua, *da;
    int *ila, *jla, *iua, *jua, *ida, *jda;  /* arrays that will be filled by
    the U, L and D in CSR format */

    if(splitAMatrix(m, problemSize, a, ia, ja, &la, &ila, &jla, &ua, &iua, &jua,
    &da, &ida, &jda)){
      printf("ERROR : Splitting matrix failed.\n");
      return EXIT_FAILURE;
    }

    if (PREC_DEBUG){
      // code is in preconditionner debugging mode for SGS solving
      printf("DEBUGGING : preconditionning for SGS : \n");
      printf("Printing the splitting of matrix A\n");
      printf("Printing the arrays...\n\n");
      printSplitSGSArrays(la, ua, ila, iua, jla, jua, da, ida, jda,
        problemSize);
    }

    // Inversing U and L

    double *invLa, *invUa;
    int *invJla, *invIla, *invJua, *invIua; /* arrays that will be filled with
    L^-1 and U^-1 in CRS format */

    int LOW = 0; // this is a lower triangular matrix
    if(inverseMatrix(problemSize, &invLa, &invJla, &invIla, la, jla, ila, LOW)){
      printf("ERROR : inversing of LA failed.\n");
      free(la); free(ua); free(da); free(ila); free(jla); free(iua); free(jua);
      free(invLa); free(invJla); free(invIla); free(ida); free(jda);
      return EXIT_FAILURE;
    }

    int UP = 1; // this is a upper triangular matrix
    if(inverseMatrix(problemSize, &invUa, &invJua, &invIua, ua, jua, iua, UP)){
      printf("ERROR : inversing of UA failed.\n");
      free(la); free(ua); free(da); free(ila); free(jla); free(iua); free(jua);
      free(invLa); free(invUa); free(invJla); free(invIla); free(invJua);
      free(invIua); free(ida); free(jda);
      return EXIT_FAILURE;
    }

    free(la); free(ua); free(jla); free(jua); free(ila); free(iua); /* freeing
    memory, only inverse matrix are needed */

    if (PREC_DEBUG){
      // code is in preconditionner debugging mode for SGS solving
      printf("DEBUGGING : inversion of matrix for SGS\n");
      printf("Printing the arrays...\n\n");
      printInvSGSArrays(invLa, invUa, invIla, invIua, invJla, invJua,
        problemSize);
    }

      // solution arrays

    double *x0 = malloc(problemSize * sizeof(double));
    *x = malloc(problemSize * sizeof(double)); /* allocation of memory
    of x_m and x_{m+1} */

    if (x0 == NULL || *x == NULL){
      printf("ERROR : not enough to initialise the solutions arrays\n");
      //free(prec); free(jPrec); free(iPrec);
      return EXIT_FAILURE;
    }

    for (int i = 0; i < problemSize; i++){
      // initialising the solution vectors to null vector
      x0[i] = 0.0;
      (*x)[i] = 0.0;
    }

      // residue

    double *residue = malloc(problemSize * sizeof(double)); /* memory allocation
    for the residue */

    if (residue == NULL){
      printf("ERROR : not enough memory to initialise the residue array\n");
      free(x0); free(*x);
      return EXIT_FAILURE;
    }

    if (getResidue(problemSize, b, a, ja, ia, x0, &residue)){
      printf("ERROR : could not compute the initial residue\n");
      free(x0); free(*x); free(residue);
      return EXIT_FAILURE;
    }

      // correction

    double *correction = malloc(problemSize * sizeof(double)); /* memory
    allocation for the corretcion array d_m */

    if (correction == NULL){
      printf("ERROR : not enough memory to create correction array\n");
      free(x0); free(*x); free(residue);
      return EXIT_FAILURE;
    }

    // Iteration loop

    int iter = 0; // number of iteration

    while(iter <= maxIterations &&
      normVector(problemSize, residue) > minResidue){
        /* loops as long as the residue is greater than minResidue and number of
        iteration is smaller than maxIterations */

        // computing correction

        if(matrixVectorMultCSR(problemSize, invLa, invJla, invIla, residue,
            &correction)){
          printf("ERROR : could not compute correction (LA^-1)\n");
          free(x0); free(*x); free(residue); free(correction);
          return EXIT_FAILURE;
        }

        if(matrixVectorMultCSR(problemSize, da, jda, ida, residue,
            &correction)){
          printf("ERROR : could not compute correction (DA)\n");
          free(x0); free(*x); free(residue); free(correction);
          return EXIT_FAILURE;
        }

        if(matrixVectorMultCSR(problemSize, invUa, invJua, invIua, residue,
            &correction)){
          printf("ERROR : could not compute correction (UA^-1)\n");
          free(x0); free(*x); free(residue); free(correction);
          return EXIT_FAILURE;
        }

        // adding the correction

        if(vectorAddition(problemSize, x0, correction, x)){
          printf("ERROR : could not add correction\n");
          free(x0); free(*x); free(residue); free(correction);
          return EXIT_FAILURE;
        }

        // computing new residue

        if (getResidue(problemSize, b, a, ja, ia, *x, &residue)){
          printf("ERROR : could not compute the residue\n");
          free(x0); free(*x); free(residue); free(correction);
          return EXIT_FAILURE;
        }

        // making x_{m+1} become x_m

        for (int i = 0; i < problemSize; i++){
          // iterating through the solutions
          x0[i] = (*x)[i];
        }

        // incrementing iteration

        iter ++;
      }

      printf("Number of iterations required to converge : %d\n", iter - 1);
      printf("Norm of last residue : %f\n\n", normVector(problemSize, residue));

      free(residue); free(x0); free(invLa); free(invUa); free(da); free(invJla);
      free(invJua); free(jda); free(invIla); free(invIua); free(ida);
      free(correction); // freeing memory

    return EXIT_SUCCESS;
  }
