// Author : Denis Verstraeten
// Created on : 23/11/2017

/* this module has been rewritten, although it is widely similar to the
provided code from 2016-2017 project */

#include "linear_problem_generator.h"

int generate_problem(int m, double L, double step, int *problemSize, \
  int **ia, int **ja, double **a, double **b, double **dirichletCond){

  /* Initialization of the parameters of the problem */

  int xNumber = m - 1; // number of points in the x direction
  int yNumber = m; // number of points in the y direction
  *problemSize = xNumber * yNumber; // number of unknowns of the problem
  int nnz = 5 * (*problemSize) - 2 * xNumber - 2 * yNumber; /* non zero element
  of the matrix A */

  /* Allocation of memory for the sparse matrix A using CSR format and
  the vector b*/

  *ia = malloc((*problemSize + 1) * sizeof(int)); /* allocation of memory for
  the ia array of the matrix A */
  *ja = malloc(nnz * sizeof(int)); /* allocation of memory for the ja array of
  the matrix A */
  *a = malloc(nnz * sizeof(double)); /* allocation of memory for the a array of
  the matrix A */
  *b = malloc(*problemSize * sizeof(double)); /* allocation of memory for the
  array of the vector b */
  *dirichletCond = malloc(*problemSize * sizeof(double)); /* allocation of
  memory for the Dirichlet condition vector */

  /* checks whether the allocation of memory was successful, returns an error
  if not */

  if(*ia == NULL || *ja ==  NULL || *a == NULL || *b == NULL || \
  *dirichletCond == NULL){
    printf("\n ERROR : not enough memory to generate the problem\n\n");
    return EXIT_FAILURE;
  }

  /* main part of this module, takes care of filling the arrays */

  nnz = 0; /* nnz is back to 0, because it will be used and incremented to fill
  the arrays */
  int equationNumber; /* this is the index of the point using the lexicographic
  numbering */

  for(int iy = 0; iy < yNumber; iy++){
    for(int ix = 0; ix < xNumber; ix++){
      /* this loop goes through every point of the problem using the
      lexicographic numbering */
      equationNumber = ix + iy * xNumber;
      (*ia)[equationNumber] = nnz; /* this indicates the start of a new line in
      the ia array */

      /* filling the independent b array by default equal to 0 */

      (*b)[equationNumber] = 0.0;

      /* filling the Dirichlet condition array */

      double y = step * iy; // y coordinate of the point
      (*dirichletCond)[equationNumber] = dirichletCondValue(L, y);

      /* starts to fill in the matrix */

      /* south neighbour element */

      if (iy > 0){
        // element not on south edge
        if (ix == 0){
          // element on west edge => Von Neumann condition : dT/dx = 0
          (*a)[nnz] = 0.5;
          (*ja)[nnz] = equationNumber - xNumber;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = 1.0;
          (*ja)[nnz] = equationNumber - xNumber;
          nnz++;
        }
      }

      /* west neighbour element */

      if (ix != 0){
        // element not on west edge
        if (iy == 0 || iy == yNumber - 1){
          // element on south or north edge => Von Neumann condition : dT/dy = 0
          (*a)[nnz] = 0.5;
          (*ja)[nnz] = equationNumber - 1;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = 1.0;
          (*ja)[nnz] = equationNumber - 1;
          nnz++;
        }
      }

      /* diagonal element */

      if (iy == 0 || iy == yNumber - 1){
        // element on south or north edge
        if (ix == 0){
          /* element on the south west or north west corner
          => Von Neumann condition : dT/dx = dT/dy = 0 */
          (*a)[nnz] = -1.0;
          (*ja)[nnz] = equationNumber;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = -2.0;
          (*ja)[nnz] = equationNumber;
          nnz++;
        }
      } else {
        // element not on south or north edge
        if (ix == 0){
          // element on west edge => Von Neumann condition : dT/dx = 0
          (*a)[nnz] = -2.0;
          (*ja)[nnz] = equationNumber;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = -4.0;
          (*ja)[nnz] = equationNumber;
          nnz++;
        }
      }

      /* east neighbour element */

      if (ix < xNumber - 1){
        // element not on the east edge
        if (iy == 0 || iy == yNumber - 1){
          // element on the north or east edge
          (*a)[nnz] = 0.5;
          (*ja)[nnz] = equationNumber + 1;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = 1.0;
          (*ja)[nnz] = equationNumber + 1;
          nnz++;
        }
      } else {
        // element on the east edge=>Dirchlet condition T = exp(sqrt(1+(y/L)^2))
        if (iy == 0 || iy == yNumber - 1){
          // element on the north or south east corner
          (*b)[equationNumber] = - 0.5 * (*dirichletCond)[equationNumber];
        } else {
          // element on the east edge
          (*b)[equationNumber] = - (*dirichletCond)[equationNumber];
        }
      }

      /* north neighbour element */

      if (iy < yNumber - 1){
        // element not the north edge
        if (ix == 0) {
          //element on the west edge
          (*a)[nnz] = 0.5;
          (*ja)[nnz] = equationNumber + xNumber;
          nnz++;
        } else {
          // internal element
          (*a)[nnz] = 1.0;
          (*ja)[nnz] = equationNumber + xNumber;
          nnz++;
        }
      }
    }
  }

  (*ia)[equationNumber + 1] = nnz; /* not really a value of ia, just a saving
  of nnz */
  printf("Linear problem successfully generated\n\n");
  return EXIT_SUCCESS; // usual function return

}
