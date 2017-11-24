// Author : Denis Verstraeten

/* this module has been rewritten, although it is widely similar to the provided
code from 2016-2017 project */

#include "linear_problem_generator.c"

int generate_problem(int m, double L, int *problemSize, int **ia, int **ja, double **a, double **b){

  /* Initialization of the parameters of the problem */

  int xNumber = m - 1; // number of points in the x direction
  int yNumber = m; // number of points in the y direction
  *problemSize = xNumber * yNumber; // number of unknowns of the problem
  double step = L / (m - 1); // length of the discretization step
  int nnz = 5 * (*problemSize) - 2 * xNumber - 2 * yNumber; // non zero elements of the matrix A

  /* Allocation of memory for the sparse matrix A using CSR format and the vector b*/

  *ia = malloc((*problemSize + 1) * sizeof(int)); // allocation of memory for the ia array of the matrix A
  *ja = malloc(nnz * sizeof(int)); // allocation of memory for the ja array of the matrix A
  *a = malloc(nnz * sizeof(double)); // allocation of memory for the a array of the matrix A
  *b = malloc(*problemSize * sizeof(double)) // allocation of memory for the array of the vector b

  /* checks whether the allocation of memory was successful, returns an error if not */

  if(*ia == NULL || *ja ==  NULL || *a == NULL){
    printf("\n ERROR : not enough memory to generate the problem\n\n");
    return 1;
  }

  /* main part of this module, takes care of filling the arrays */

  nnz = 0; // nnz is back to 0, because it will be used and incremented to fill the arrays
  int equationNumber; // this is the index of the point using the lexicographic numbering

  for(int iy = 0; iy < yNumber, iy++){
    for(int ix = 0; ix < xNumber, ix++){
      /* this loop goes through every point of the problem using the lexicographic numbering */
      equationNumber = ix + iy * xNumber;
      (*ia)[nnz] = equationNumber; // this indicates the start of a new line in the ia array

      /* filling the independent b array by default equal to 0 */

      (*b)[equationNumber] = 0.0;

      /* starts to fill in the matrix */

      /* south neighbour elements */

      if (iy > 0 && ix == 0){
        // elements on the upper or lower edge => Von Neumann condition : dT/dy = 0
        (*a)[nnz] = 0.5;
        (*ja)[nnz] = equationNumber - xNumber;
        nnz++;
      } else if (iy > 0 && ix > 0){
        // internal elements
        (*a)[nnz] = 1.0;
        (*ja)[nnz] = equationNumber - xNumber;
        nnz++;
      }

      /* west neighbour elements */

      if ((ix > 0 && iy == 0) || (ix > 0 && iy == yNumber - 1)){
        // elements on the upper or lower edge => Von Neumann condition : dT/dy = 0
        (*a)[nnz] = 0.5;
        (*ja)[nnz] = equationNumber - 1;
        nnz++;
      } else if (ix > 0 && (iy > 0 && iy < yNumber - 1)){
        // internal elements
        (*a)[nnz] = 1.0;
        (*ja)[nnz] = equationNumber - 1;
        nnz++;
      }

      /* diagonal elements */

      if ((ix == 0 && iy == 0) || (ix == 0 && iy == yNumber - 1)){
        // elements on the lower or upper west corner => Von Neumann condition : dT/dx = dT/dy = 0
        (*a)[nnz] = -1.0;
        (*ja)[nnz] = equationNumber;
        nnz++;
      } else if ((iy == 0 && ix > 0) || (iy == yNumber - 1 && ix > 0)){
        // elements on the upper or lower edge => Von Neumann condition : dT/dy = 0
        (*a)[nnz] = -2.0;
        (*ja)[nnz] = equationNumber;
        nnz++;
      } else if (ix == 0 && (iy > 0 && iy < yNumber - 1)){
        // elements on the west edge => Von Neumann condition : dT/dx = 0
        (*a)[nnz] = -2.0;
        (*ja)[nnz] = equationNumber;
        nnz++;
      } else {
        // internal elements
        (*a)[nnz] = -4.0;
        (*ja)[nnz] = equationNumber;
        nnz;
      }

      /* east neighbour elements */



    }
  }
  printf("Linear problem successfully generated\n");
  return 0;
}
