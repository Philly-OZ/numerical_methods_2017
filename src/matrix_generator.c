// Author : Denis Verstraeten
// Created on : 23/11/2017

/* this module has been rewritten, although it is widely similar to the
provided code from 2016-2017 project */

#include "matrix_generator.h"

int generateProblem(int m, double L, double step, int *problemSize,
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

  if(*ia == NULL || *ja ==  NULL || *a == NULL || *b == NULL ||
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

  (*ia)[equationNumber + 1] = nnz; // saving of nnz
  printf("Linear problem successfully generated\n");
  return EXIT_SUCCESS; // usual function return

}

int generateProlongation(int mFineGrid, int mCoarseGrid,
  double **prol, int **jProl, int **iProl){

    /* this function returns the arrays corresponding to the prolongation
    matrix P that will be used for the multi-grid scheme.
    The matrix P is of dimension nxm with n > m */

    /* Parameters */

    // Generating the fine grid parameters

    int fineYNumber = mFineGrid; /* number of points in the y direction relative
    to the fine grid */
    int fineXNumber = mFineGrid - 1; /* number of points in the x direction
    relative to the fine grid */
    int fineProblemSize = fineXNumber * fineYNumber; /* number of equations and
    unknown in the fine problem */

    // Generating the coarse grid parameters

    int coarseYNumber = mCoarseGrid; /* number of points in the y direction
    relative to the coarse grid */
    int coarseXNumber = mCoarseGrid - 1; /* number of points in the x direction
    relative to the coarse grid */
    int coarseProblemSize = coarseXNumber * coarseYNumber; /* number of equations
    and unknown in the coarse problem */

    /* Memory allocation */

    // Computing non zero elements of P

    int nnz = coarseProblemSize + coarseYNumber + 2 * coarseYNumber *
    (coarseXNumber - 1) + 2 * (coarseXNumber + 1) * (coarseYNumber - 1) +
    4 * (coarseXNumber - 1) * (coarseYNumber - 1); /* number of non zero
    elements in the prolongation matrix */

    // Allocating memory

    *prol = malloc(nnz * sizeof(double));
    *jProl = malloc(nnz * sizeof(int));
    *iProl = malloc((fineProblemSize + 1) * sizeof(int)); /* allocation of
    memory for the arrays describing P */

    // Allocation check

    if (prol == NULL){
      printf("ERROR : not enough memory to generate prol array\n");
      return EXIT_FAILURE;
    }

    if(jProl == NULL){
      printf("ERROR : not enough memory to generate jProl array\n");
      free(prol);
      return EXIT_FAILURE;
    }

    if (iProl == NULL){
      printf("ERROR : not enough memory to generate iProl array\n");
      free(prol); free(jProl);
      return EXIT_FAILURE;
    }

    /* Filling in the arrays */

    nnz = 0; // back to zero, it will be used and incremented to fill the arrays
    int equationNumber; /* this is the index of the point using the
    lexicographic numbering */

    for (int i = 0; i < fineYNumber; i++){
      // iterating through the y direction in the fine grid
      for (int j = 0; j < fineXNumber; j++){
        // iterating through the x direction in the fine grid

        equationNumber = i * fineXNumber + j; // number of the current equation
        (*iProl)[equationNumber] = nnz;

        if (i % 2 == 0){
          // the point is on a line made of coarse points
          if(j < fineXNumber - 1){
            // the point is not on the East border of the fine grid
            if (j % 2 == 0){
              // the point is on a coarse point
              (*prol)[nnz] = 1;
              (*jProl)[nnz] = (i / 2) * coarseXNumber + (j / 2);
              nnz ++;
            } else {
              // the point is between two coarse points in the x direction
              (*prol)[nnz] = 0.5;
              (*jProl)[nnz] = (i / 2) * coarseXNumber + (int) j / 2;
              nnz ++;
              (*prol)[nnz] = 0.5;
              (*jProl)[nnz] = (i / 2) * coarseXNumber + (int) ((j / 2) + 1);
              nnz ++;
            }
          } else {
            // the point is on the on the East border of the fine grid
            (*prol)[nnz] = 0.5;
            (*jProl)[nnz] = (i / 2) * coarseXNumber + (int) j / 2;
            nnz ++;
          }
        } else {
          // the point is not on a line made of coarse points
          if(j < fineXNumber - 1){
            // the point is not on the East border of the fine grid
            if (j % 2 == 0){
              // the point is between two coarse points in the y direction
              (*prol)[nnz] = 0.5;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (j / 2);
              nnz ++;
              (*prol)[nnz] = 0.5;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (j / 2) +
                coarseXNumber;
              nnz ++;
            } else {
              // the point is in the middle on four coarse points
              (*prol)[nnz] = 0.25;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2);
              nnz ++;
              (*prol)[nnz] = 0.25;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2)
              + 1;
              nnz ++;
              (*prol)[nnz] = 0.25;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2)
              + coarseXNumber;
              nnz ++;
              (*prol)[nnz] = 0.25;
              (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2)
              + coarseXNumber + 1;
              nnz ++;
            }
          } else {
            /* the point is on the East border of the fine grid, between two
            coarse points NW and SW */
            (*prol)[nnz] = 0.25;
            (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2);
            nnz ++;
            (*prol)[nnz] = 0.25;
            (*jProl)[nnz] = ((int) (i / 2) * coarseXNumber) + (int) (j / 2)
            + coarseXNumber;
            nnz ++;
          }
        }
      }
    }

    (*iProl)[equationNumber + 1] = nnz; // saving of nnz
    return EXIT_SUCCESS;

  }
