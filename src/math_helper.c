// Author : Denis Verstraeten
// Created on : 27/11/2017

/* this module is used for all calculations */

#include "math_helper.h"

int factorisation(int problemSize, double *a, int *ja, int *ia, void **Numeric){
  // this factorises a matrix to solve it
	int statut;
	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
	void *Symbolic;
	statut = umfpack_di_symbolic(problemSize, problemSize, ia, ja, a, &Symbolic,
     Control, Info);
	if (statut < 0)
	{
		printf("\nERROR : UMF_PACK_SYMBOLIC FAILED");
		return EXIT_FAILURE;
	}
	statut = umfpack_di_numeric(ia, ja, a, Symbolic, Numeric, Control, Info);
	if (statut < 0)
	{
		printf("ERROR : UMF_PACK_NUMERIC FAILED\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

int resolution(int *ia, int *ja, double *a, double *x, double *b,
  void **Numeric){
  // this will solve the linear system Ax=b with the matrix A already factorised
	int statut;
	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
	statut = umfpack_di_solve(UMFPACK_At, ia, ja, a, x, b, *Numeric, Control,
     Info);
	if (statut < 0)
	{
		printf("ERROR : UMF_PACK_SOLVE FAILED\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

double dirichletCondValue(double L, double y){
	/* this function returns the value of the Dirichlet boundary condition at
	point y */
	return exp(sqrt(1+square(y / L)));
}

int splitAMatrix(int m, int problemSize, double *a, int *ia, int *ja,
	double **la, int **ila, int **jla, double **ua, int **iua, int **jua,
	double **da){
	/* this funtions breaks the A matrix into LA, DA, UA so that
	A = LA + UA - DA.
	LA :  lower triangular part of A
	UA : upper triangular part of A
	DA : diagonal elements of A */

	/* Initialization of the parameters of the problem */

	int xNumber = m - 1; // number of points in the x direction
  int yNumber = m; // number of points in the y direction
	int nnzLA = 3 * problemSize - xNumber - yNumber; /* number of non zero
	elements of the la array corresponding to LA matrix */
	int nnzUA = 3 * problemSize - xNumber - yNumber; /* number of non zero
	elements of the la array corresponding to UA matrix */

	/* Allocation of memory for the CSR arrays corresponding to LA and UA matrix
	and the array containing DA */

	// LA matrix

	*la = malloc(nnzLA * sizeof(double)); /* allocation of memory for the la array
	corresponding to the LA matrix */
	*jla = malloc(nnzLA * sizeof(int)); /* allocation of memory for the jla array
	corresponding to the LA matrix */
	*ila = malloc((problemSize + 1) * sizeof(int)); /* allocation of memory of the
	ila array corresponding to the LA matrix */

	// UA matrix

	*ua = malloc(nnzUA * sizeof(double)); /* allocation of memory for the ua array
	corresponding to the UA matrix */
	*jua = malloc(nnzUA * sizeof(int)); /* allocation of memory for the jua array
	corresponding to the UA matrix */
	*iua = malloc((problemSize + 1) * sizeof(int)); /* allocation of memory of the
	iua array corresponding to the UA matrix */

	// DA matrix : this matrix is only stored as a vector to save memory

	*da = malloc(problemSize * sizeof(double)); /* allocation of memory for the
	da array corresponding to DA matrix */

	// checks whether the allocation was successful, returns an error if not

	if (*la == NULL || *jla == NULL || *ila == NULL || *ua == NULL ||
		*jua == NULL || *iua == NULL ){
			printf("\n ERROR : not enough memory to generate preconditionning\n\n");
			return EXIT_FAILURE;
		}

	// filling in the arrays

	nnzLA = 0;
	nnzUA = 0; /* the nnz are back to 0, because they will be used and
		incremented to fill the arrays */

	int i; // making loop increment a global variable
	for (i = 0; i < problemSize; i++){
		// iterating through the lines of matrix A
		(*ila)[i] = nnzLA;
		(*iua)[i] = nnzUA;
		for (int j = ia[i]; j < ia[i + 1]; j++){
			// iterating through the non zero elements of the ith line
			if (ja[j] < i){
				// strictly lower triangular part of matrix A
				(*la)[nnzLA] = a[j];
				(*jla)[nnzLA] = ja[j];
				nnzLA ++;
			} else if (ja[j] == i){
				// diagonal part of matrix A
				(*la)[nnzLA] = a[j];
				(*jla)[nnzLA] = ja[j];
				nnzLA ++;
				(*ua)[nnzUA] = a[j];
				(*jua)[nnzUA] = ja[j];
				nnzUA ++;
				(*da)[i] = a[j];
			} else if (ja[j] > i){
				// strictly upper triangular part of matrix A
				(*ua)[nnzUA] = a[j];
				(*jua)[nnzUA] = ja[j];
				nnzUA ++;
			}
		}
	}
	(*ila)[i] = nnzLA; // saving of nnz of LA
	(*iua)[i] = nnzUA; // saving of nnz of UA
	return EXIT_SUCCESS;
}
