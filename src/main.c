// Author : Denis Verstraeten
// Created on : 24/11/2017

/* this is the main module of the project */

#include "main.h"

int main(int argc, char **argv){

  if (argc <= 5){

    // correct number of arguments

    // option variables declaration

    int DEBUG_A_MATRIX; /* boolean variable indicating whether the code is in
    A matrix debugging mode */
    int DEBUG_LINEAR_SYSTEM; /* boolean variable indicating whether the code is
    in linear system debugging mode */
    int UMF_SOLVE; /* boolean variable indicating whether the code is solving
    the linear problem using UMF Pack */
    int PRINT_SOLUTION; /* boolean variable indicating whether the code needs
    to print the solution of the solved problem */
    int PLOT_SOLUTION; /* boolean variable indicating whether the code needs to
    plot the solution */
    generateOptions(&DEBUG_A_MATRIX, &DEBUG_LINEAR_SYSTEM, &UMF_SOLVE, \
      &PRINT_SOLUTION, &PLOT_SOLUTION, argc, argv);

    // variables declaration

    int m = 1000; // number of points in the y direction
    double L = 0.2; // size of the square membrane
    double step = L / (m - 1); // length of the discretization step
    int problemSize, *ia, *ja; /* number of unknowns of the problem, arrays that
     will be filled by the matrix in CSR format */
    double *a, *b; /* array containing the non zero elements of the matrix of
    the problem in CSR format, array of independent terms */
    double *dirichletCond; /* array containing the Dirichlet condition of the
    East edge */

    // Program starting

    printf("\n~~~~~~~~~~~~~~~~~~~\n");
    printf("PROGRAM STARTING...\n");
    printf("~~~~~~~~~~~~~~~~~~~\n\n");

    if (generate_problem(m, L, step, &problemSize, &ia, &ja, &a, &b, \
      &dirichletCond)){
      /* this will end the program if there was a problem with the creation of
      the arrays */
      return EXIT_FAILURE;
    }

    if (DEBUG_A_MATRIX) {
      // if debug a matrix is enabled
      printf("A matrix debugging is enabled\n");
      printAArrays(a, ja, ia, problemSize);
    }

    if (DEBUG_LINEAR_SYSTEM){
      // if debug linear system is enabled
      printf("Linear system debugging is enabled\n");
      printLinearSystemArrays(a, ja, ia, b, problemSize);
    }

    if (UMF_SOLVE){
      // if UMF Pack solving of the linear system is enabled
      printf("Solving the problem using UMF Pack...\n\n");
      double *T = malloc(problemSize * sizeof(double)); /* vector containing the
      solution computed using UMF Pack */
      if (umfSolve(problemSize, a, ja, ia, T, b)){
        printf("ERROR : UMF Pack solving failed\n");
        return EXIT_FAILURE;
      }
      if (PRINT_SOLUTION){
        printf("Printing the solution obtained with UMF Pack\n\n");
        // if print solution mode is enabled
        for (int i = 0; i < problemSize; i++){
          // iterating through the solution vector
          printf("%f\n", T[i]);
        }
        printf("\n");
      }
      if (PLOT_SOLUTION){
        printf("Saving the plot of the solution obtained with UMF Pack\n");
        // if plot solution mode is enabled
        if (plot(m, step, T, dirichletCond)) {
          printf("ERROR : could not plot the solution\n");
          return EXIT_FAILURE;
        }
        printf("The result is saved in graphics directory\n\n");
      }
    }

    printf("Program ending...\n\n");
    return EXIT_SUCCESS;
  } else {
    // too many arguments
    printf("Too many arguments passed, program quitting...\n");
    return EXIT_FAILURE;
  }
}
