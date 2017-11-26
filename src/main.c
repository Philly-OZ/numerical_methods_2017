// Author : Denis Verstraeten
// Created on : 24/11/2017

/* this is the main module of the project */

#include "main.h"

int main(int argc, char *argv[]){

  if (argc < 3){

    // correct number of arguments

    // variables declarations

    int DEBUG_A_MATRIX = argc == 2 && strcmp(argv[1], "debug=a_matrix") == 0; /* boolean
    indicating whether the code is in debug=a_matrix mode */
    int m = 4; // number of points in the y direction
    double L = 0.2; // size of the square membrane
    int problemSize, *ia, *ja; /* number of unknowns of the problem, arrays that will be
    filled by the matrix in CSR format */
    double *a, *b; /* array contzining the non zero elements of the matrix of the
    problem in CSR format, array of independent terms */

    printf("\n~~~~~~~~~~~~~~~~~~~\n");
    printf("PROGRAM STARTING...\n");
    printf("~~~~~~~~~~~~~~~~~~~\n\n");

    if (generate_problem(m, L, &problemSize, &ia, &ja, &a, &b)){
      /* this will end the program if there was a problem with the creation of the
      arrays */
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  } else {
    // too many arguments
    printf("Too many arguments passed, program quitting...\n");
    return EXIT_FAILURE;
  }
}
