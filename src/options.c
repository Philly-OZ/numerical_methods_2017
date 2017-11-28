// Author : Denis Verstraeten
// Created on : 26/11/2017

/* This module contains all the functions that are called during a debug
execution of the code */

#include "options.h"

void generateOptions(int *DEBUG_A_MATRIX, int *DEBUG_LINEAR_SYSTEM, \
  int *UMF_SOLVE, int *PRINT_SOLUTION, int *PLOT_SOLUTION, int argc, \
  char **argv){
  // generates the boolean values corresponding to the selected options
  *DEBUG_A_MATRIX = 0;
  *DEBUG_LINEAR_SYSTEM = 0;
  *UMF_SOLVE = 0;
  *PRINT_SOLUTION = 0;
  *PLOT_SOLUTION = 0;
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
    }
  }
}

void printAArrays(double *a, int *ja, int *ia, int problemSize){
  // prints the values contained the arrays corresponding to A matrix
  printf("Printing the CSR arrays corresponding to the A matrix\n\n");
  printf("a array :\n");
  printf("---------\n");
  for (int i = 0; i < ia[problemSize]; i++){
    // iterating through elements of the a array
    printf("%f\n", a[i]);
  }
  printf("\n\nja array :\n");
  printf("----------\n");
  for (int i = 0; i < ia[problemSize]; i++){
    // iterating through elements of the ja array
    printf("%i\n", ja[i]);
  }
  printf("\n\nia array :\n");
  printf("----------\n");
  for (int i = 0; i < problemSize + 1; i++){
    // iterating through elements of the ia array
    printf("%i\n", ia[i]);
  }
  printf("\n\n------------\n\n");
}

void printLinearSystemArrays(double *a, int *ja, int *ia, double *b, \
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
