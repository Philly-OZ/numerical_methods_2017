// Author : Denis Verstraeten
// Created on : 24/11/2017

/* this is the main module of the project */

#include "main.h"

int main(int argc, char **argv){

  if (argc <= 3){

    // correct number of arguments

    // debugging variables declaration

    int DEBUG_A_MATRIX; /* boolean variable indicating whether the code is in
    A matrix debugging mode */
    int DEBUG_LINEAR_SYSTEM; /* boolean variable indicating whether the code is
    in linear system debugging mode */
    generateOptions(&DEBUG_A_MATRIX, &DEBUG_LINEAR_SYSTEM, argc, argv);

    // variables declaration

    int m = 4; // number of points in the y direction
    double L = 0.2; // size of the square membrane
    int problemSize, *ia, *ja; /* number of unknowns of the problem, arrays that will be
    filled by the matrix in CSR format */
    double *a, *b; /* array containing the non zero elements of the matrix of the
    problem in CSR format, array of independent terms */

    // Program starting

    printf("\n~~~~~~~~~~~~~~~~~~~\n");
    printf("PROGRAM STARTING...\n");
    printf("~~~~~~~~~~~~~~~~~~~\n\n");

    if (generate_problem(m, L, &problemSize, &ia, &ja, &a, &b)){
      /* this will end the program if there was a problem with the creation of the
      arrays */
      return EXIT_FAILURE;
    }

    if (DEBUG_A_MATRIX) {
      // if debug a matrix is enabled
      printf("A matrix debugging is enabled\n");
      printAArrays(a, ja, ia, problemSize);
    }

    if(DEBUG_LINEAR_SYSTEM){
      // if debug linear system is enabled
      printf("Linear system debugging is enabled\n");
      printLinearSystemArrays(a, ja, ia, b, problemSize);
    }

    printf("Program ending...\n\n");
    return EXIT_SUCCESS;
  } else {
    // too many arguments
    printf("Too many arguments passed, program quitting...\n");
    return EXIT_FAILURE;
  }
}
