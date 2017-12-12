// Author : Denis Verstraeten
// Created on : 26/11/2017

/* This module contains all the functions that are called during a debug
execution of the code */

#include "options.h"

void generateOptions(int *DEBUG_A_MATRIX, int *DEBUG_LINEAR_SYSTEM,
  int *UMF_SOLVE, int *PRINT_SOLUTION, int *PLOT_SOLUTION, int *SGS_SOLVE,
  int *PREC_DEBUG, int argc, char **argv){
  // generates the boolean values corresponding to the selected options
  *DEBUG_A_MATRIX = 0;
  *DEBUG_LINEAR_SYSTEM = 0;
  *UMF_SOLVE = 0;
  *PRINT_SOLUTION = 0;
  *PLOT_SOLUTION = 0;
  *SGS_SOLVE = 0;
  *PREC_DEBUG = 0;
  for (int i = 1; i < argc; i++){
    // iterating through every passed arguments
    if (strcmp(argv[i], "-a") == 0){
       // if in debugging A matrix mode
      *DEBUG_A_MATRIX = 1;
    } else if (strcmp(argv[i], "-linsys") == 0){
      // if in linear system debugging mode
      *DEBUG_LINEAR_SYSTEM = 1;
    } else if (strcmp(argv[i], "-umf-solve") == 0){
      // if the system is solved using UMF Pack
      *UMF_SOLVE = 1;
    } else if (strcmp(argv[i], "-print-s") == 0){
      // if the solution needs to be printed
      *PRINT_SOLUTION = 1;
    } else if (strcmp(argv[i], "-plot") == 0){
      // if the solution needs to be plotted
      *PLOT_SOLUTION = 1;
    } else if (strcmp(argv[i], "-sgs-solve") == 0){
      // if the system is solved using SGS
      *SGS_SOLVE = 1;
    } else if (strcmp(argv[i], "-sgs-prec-debug") == 0){
      // if in debugging preconditionning for SGS mode
      *PREC_DEBUG = 1;
    }
  }
}

void printAArrays(double *a, int *ja, int *ia, int problemSize){
  // prints the values contained the arrays corresponding to A matrix
  printf("Printing the CSR arrays corresponding to the A matrix\n\n");
  printf("a and ja arrays :\n");
  printf("-----------------\n");
  for (int i = 0; i < ia[problemSize]; i++){
    // iterating through elements of the a and ja arrays
    printf("a[%d] = %f, ja[%d] = %d\n", i, a[i], i, ja[i]);
  }
  printf("\n\nia array :\n");
  printf("----------\n");
  for (int i = 0; i < problemSize + 1; i++){
    // iterating through elements of the ia array
    printf("ia[%d] = %i\n", i, ia[i]);
  }
  printf("\n\n------------\n\n");

  printf("Number of nnz elements : \n");
  printf("A : %d\n", ia[problemSize]);

  printf("\n\n------------\n\n");
}

void printLinearSystemArrays(double *a, int *ja, int *ia, double *b,
  int problemSize){
  /* prints the values contained the arrays corresponding to A matrix and the
  independent term b array */
  printAArrays(a, ja, ia, problemSize);
  printf("Printing the array corresponding to b vector\n\n");
  printf("b array :\n");
  printf("---------\n");
  for (int i = 0; i < problemSize; i++){
    // iterating through b array
    printf("%f\n", b[i]);
  }
  printf("\n\n------------\n\n");

}

void printSplitSGSArrays(double *la, double *ua, int *ila, int *iua, int *jla,
  int *jua, double *da, int *ida, int *jda, int problemSize){
    /* prints the values contained in the arrays representing L, U and D
    matrixes used in SGS preconditionning */
    int nnzLA = ila[problemSize]; // number of nnz in the L matrix
    int nnzUA = iua[problemSize]; // numerb of nnz in the U matrix
    printf("Printing the arrays corresponding to matrix L \n\n");
    printf("la and jla arrays : \n");
    printf("--------------------\n\n");
    for (int i = 0; i < nnzLA; i++){
      // iterating though elements of la and jla arrays
      printf("la[%d] = %f, jla[%d] = %d\n", i, la[i], i, jla[i]);
    }
    printf("\n\nila array : \n");
    printf("------------\n\n");
    for (int i = 0; i < problemSize + 1; i++){
      // iterating though elements of ila array
      printf("ila[%d] = %d\n", i, ila[i]);
    }
    printf("\n\n------------\n\n");

    printf("Printing the arrays corresponding to matrix U \n\n");
    printf("ua and jua arrays : \n");
    printf("--------------------\n\n");
    for (int i = 0; i < nnzUA; i++){
      // iterating though elements of la and jla arrays
      printf("ua[%d] = %f, jua[%d] = %d\n", i, ua[i], i, jua[i]);
    }
    printf("\n\niua array : \n");
    printf("------------\n\n");
    for (int i = 0; i < problemSize + 1; i++){
      // iterating though elements of ila array
      printf("iua[%d] = %d\n", i, iua[i]);
    }
    printf("\n\n------------\n\n");

    printf("Printing the array corresponding to matrix D \n\n");
    printf("da and jda arrays : \n");
    printf("-------------------\n\n");
    for (int i = 0; i < problemSize; i++){
      // iterating through da array
      printf("da[%d] = %f, jda[%d] = %d\n", i, da[i], i, jda[i]);
    }
    printf("\n\nida array : \n");
    printf("------------\n\n");
    for (int i = 0; i < problemSize + 1; i++){
      // iterating though elements of ida array
      printf("ida[%d] = %d\n", i, ida[i]);
    }
    printf("\n\n------------\n\n");

    printf("Number of nnz elements : \n");
    printf("LA : %d, UA : %d, DA : %d\n\n", ila[problemSize], iua[problemSize],
      ida[problemSize]);

  }

void printInvSGSArrays(double *invLa, double *invUa, int *invIla, int *invIua,
   int *invJla, int *invJua, int problemSize){
    /* prints the values contained in the arrays representing L^-1 and U^-1
    matrixes used in SGS preconditionning */
    int nnzInvLA = invIla[problemSize]; // number of nnz in the L matrix
    int nnzInvUA = invIua[problemSize]; // numerb of nnz in the U matrix
    printf("Printing the arrays corresponding to matrix L^-1 \n\n");
    printf("invLa and invJla arrays : \n");
    printf("--------------------\n\n");
    for (int i = 0; i < nnzInvLA; i++){
      // iterating though elements of la and jla arrays
      printf("invLa[%d] = %lf, invJla[%d] = %d\n", i, invLa[i], i, invJla[i]);
    }
    printf("\n\ninvIla array : \n");
    printf("------------\n\n");
    for (int i = 0; i < problemSize + 1; i++){
      // iterating though elements of ila array
      printf("invIla[%d] = %d\n", i, invIla[i]);
    }
    printf("\n\n------------\n\n");

    printf("Printing the arrays corresponding to matrix U^-1 \n\n");
    printf("invUa and invJua arrays : \n");
    printf("--------------------\n\n");
    for (int i = 0; i < nnzInvUA; i++){
      // iterating though elements of la and jla arrays
      printf("invUa[%d] = %lf, invJua[%d] = %d\n", i, invUa[i], i, invJua[i]);
    }
    printf("\n\ninvIua array : \n");
    printf("------------\n\n");
    for (int i = 0; i < problemSize + 1; i++){
      // iterating though elements of ila array
      printf("invIua[%d] = %d\n", i, invIua[i]);
    }
    printf("\n\n------------\n\n");

    printf("Number of nnz elements : \n");
    printf("LA^-1 : %d, UA^-1 : %d\n\n", invIla[problemSize],
      invIua[problemSize]);

    printf("\n\n------------\n\n");
  }


void printPrecSGSArrays(double *prec, int *jPrec, int *iPrec, int problemSize){
  /* prints the values contained in the arrays representing the
  preconditionner for SGS solve */
  printf("Printing the CSR arrays corresponding to the preconditionning "
    "matrix\n\n");
  printf("prec and jPrec arrays :\n");
  printf("-----------------------\n");
  for (int i = 0; i < iPrec[problemSize]; i++){
    // iterating through elements of the a and ja arrays
    printf("prec[%d] = %f, jPrec[%d] = %d\n", i, prec[i], i, jPrec[i]);
  }
  printf("\n\niPrec array :\n");
  printf("-------------\n");
  for (int i = 0; i < problemSize + 1; i++){
    // iterating through elements of the ia array
    printf("iPrec[%d] = %i\n", i, iPrec[i]);
  }
  printf("\n\n------------\n\n");

  printf("Number of nnz elements : \n");
  printf("B^-1 : %d\n", iPrec[problemSize]);

  printf("\n\n------------\n\n");

}
