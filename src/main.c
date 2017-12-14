// Author : Denis Verstraeten
// Created on : 24/11/2017

/* this is the main module of the project */

#include "main.h"

int main(int argc, char **argv){

  if (argc <= 9){

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
    int SGS_SOLVE; /* boolean variable indicating whether the code is solving
    the linear problem using the symmetric Gauss-Seidel algorithm */
    int PREC_DEBUG; /* boolean variable indicating whether the code is in debug
    mode for the preconditionning for the SGS solving */
    int NORM_RESIDUE; /* booelan variable indicating whether the code needs to
    plot the residue as a function of the number of iteration for iterative
    methods */
    generateOptions(&DEBUG_A_MATRIX, &DEBUG_LINEAR_SYSTEM, &UMF_SOLVE,
      &PRINT_SOLUTION, &PLOT_SOLUTION, &SGS_SOLVE, &PREC_DEBUG, &NORM_RESIDUE,
      argc, argv);

    // Program starting

    printf("\n~~~~~~~~~~~~~~~~~~~\n");
    printf("PROGRAM STARTING...\n");
    printf("~~~~~~~~~~~~~~~~~~~\n\n");

    // parametric number of discretised points

    int m; // number of points in the y direction, input from the user
    printf("Enter m (number of points in the y direction) : ");
    scanf("%d", &m); // saving the value of m  in the variable
    printf("This input will be modified to the closest bigger suitable"
      " value\n\n");
    int maxDepth = makeSuitableM(&m); /* corrects m to a value suitable for
    multi-grid solving and stores the maximum depth */
    if (maxDepth == EXIT_FAILURE){
      // suitable m generation failed
      return EXIT_FAILURE;
    }

    // variables declaration

    double L = 0.2; // size of the square membrane
    double step = L / (m - 1); // length of the discretization step
    int problemSize, *ia, *ja; /* number of unknowns of the problem, arrays that
     will be filled by the matrix in CSR format */
    double *a, *b; /* array containing the non zero elements of the matrix of
    the problem in CSR format, array of independent terms */
    double *dirichletCond; /* array containing the Dirichlet condition of the
    East edge */

    // generating linear problem

    printf("Generating linear problem...\n\n");

    double timeBeforeProblem = timer(); // time before generating the problem
    if (generateProblem(m, L, step, &problemSize, &ia, &ja, &a, &b,
      &dirichletCond)){
      /* this will end the program if there was a problem with the creation of
      the arrays */
      return EXIT_FAILURE;
    }
    double timeAfterProblem = timer(); // time when the problem is generated
    printf("Time taken to generate the problem : %f seconds\n\n",
  timeAfterProblem - timeBeforeProblem);
    printf("Problem features : \n");
    printf("------------------\n\n");
    printf("number of points in x direction : %d\n", m - 1);
    printf("number of points in y direction : %d (=m)\n", m);
    printf("number of unknowns (and equations) : %d\n", problemSize);
    printf("number of non zero elements of A : %d\n", ia[problemSize]);
    printf("step length : %f meters\n", step);
    printf("--------------------------\n\n");

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
      time_t timeBeforeUmfSolve = time(0);
      double timeBeforeUmfSolveMicro = timer();
      if (umfSolve(problemSize, a, ja, ia, T, b)){
        printf("ERROR : UMF Pack solving failed\n");
        free(T); free(a); free(ja); free(ia); free(b); free(dirichletCond);
        /* Realeasing memory */
        return EXIT_FAILURE;
      }
      time_t timeAfterUmfSolve = time(0);
      double timeAfterUmfSolveMicro = timer();
      double elapsedTime = timeAfterUmfSolveMicro - timeBeforeUmfSolveMicro;
      if (elapsedTime > 1){
        elapsedTime = (double) timeAfterUmfSolve - timeBeforeUmfSolve;
      }
      printf("Time taken to solve the problem using UMF Pack : %f seconds\n",
        elapsedTime );
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
        printf("Saving the plot of the solution obtained with UMF Pack...\n");
        // if plot solution mode is enabled
        if (plot(m, step, T, dirichletCond)) {
          printf("ERROR : could not plot the solution\n");
          free(T); free(a); free(ja); free(ia); free(b); free(dirichletCond);
          /* Realeasing memory */
          return EXIT_FAILURE;
        }
        printf("The result is saved in graphics directory\n\n");
      }
      free(T); // Realeasing memory
    }

    if (SGS_SOLVE){
      // if the Symmetric Gauss Seidel solving mode is enabled
      printf("Solving the problem using Symmetric Gauss Seidel...\n\n");
      double *T; // memory allocation for solution vector
      time_t timeBeforeSGSSolve = time(0);
      double timeBeforeSGSSolveMicro = timer();
      if (sgsSolve(a, ia, ja, &T, b, 1e-3, 1e4, m, problemSize, PREC_DEBUG,
          NORM_RESIDUE)){
        printf("ERROR : SGS failed\n");
        free(T); free(a); free(ja); free(ia); free(b); free(dirichletCond);
        /* Realeasing memory */
        return EXIT_FAILURE;
      }
      time_t timeAfterSGSSolve = time(0);
      double timeAfterSGSSolveMicro = timer();
      double elapsedTime = timeAfterSGSSolveMicro - timeBeforeSGSSolveMicro;
      if(elapsedTime > 1){
        elapsedTime = (double) timeAfterSGSSolve - timeBeforeSGSSolve;
      }
      printf("Time taken to solve the problem using SGS iterative"
      " method : %f seconds\n\n", elapsedTime);
      if (PRINT_SOLUTION){
        printf("Printing the solution obtained with SGS iterative"
        " method\n\n");
        // if print solution mode is enabled
        for (int i = 0; i < problemSize; i++){
          // iterating through the solution vector
          printf("%f\n", T[i]);
        }
        printf("\n");
      }
      if (PLOT_SOLUTION){
        printf("Saving the plot of the solution obtained with SGS iterative"
        " method...\n");
        // if plot solution mode is enabled
        if (plot(m, step, T, dirichletCond)) {
          printf("ERROR : could not plot the solution\n");
          free(T); free(a); free(ja); free(ia); free(b); free(dirichletCond);
          /* Realeasing memory */
          return EXIT_FAILURE;
        }
        printf("The result is saved in graphics directory\n\n");
      }
      if (NORM_RESIDUE){
        printf("Saving the plot of the residue obtained with SGS iterative"
        " method...\n");
        plotResidue("SGS iterative method", "SGS");
        printf("The result is saved in graphics directory\n\n");
        system("rm residueSGS.txt"); // removing the temporary data file
        system("mv residueSGS.png ./graphics"); /* storing the png in the
        graphics folder */
      }
      free(T); // Realeasing memory

    }

    printf("Program ending...\n\n");

    /* Releasing memory before quitting */

    free(ia); free(ja); free(a); free(b); free(dirichletCond);

    return EXIT_SUCCESS;
  } else {
    // too many arguments
    printf("Too many arguments passed, program quitting...\n");
    return EXIT_FAILURE;
  }
}
