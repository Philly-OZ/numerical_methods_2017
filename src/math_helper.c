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
	if (statut < 0){
		printf("\nERROR : UMF_PACK_SYMBOLIC FAILED");
		return EXIT_FAILURE;
	}
	statut = umfpack_di_numeric(ia, ja, a, Symbolic, Numeric, Control, Info);
	if (statut < 0){
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
	if (statut < 0){
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
			free(la); free(jla); free(ila); free(ua); free(jua); free(iua); free(da);
			// freeing memory
			return EXIT_FAILURE;
		}

	// filling in the arrays

	nnzLA = 0;
	nnzUA = 0; /* the nnz are back to 0, because they will be used and
		incremented to fill the arrays */

	for (int i = 0; i < problemSize; i++){
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
	(*ila)[problemSize] = nnzLA; // saving of nnz of LA
	(*iua)[problemSize] = nnzUA; // saving of nnz of UA
	return EXIT_SUCCESS;
}

int inverseMatrix(int problemSize, double **invA, int **invJa, int **invIa,
	double *a, int *ja, int *ia){
		/* this functions returns the inverse of the matrix A stored in CSR in
		a, ja and ia as CSR arrays invA, invJa, invIa. It uses the UMF Pack solver
		since this will actually only deals with triangular matrix.
		It solves the system Ax=b <=> x=A\b with b a vector of the canonic basis
		=> solving the system b_i gives the ith column of A^-1 */
		double *tempInvA = malloc(square(problemSize) * sizeof(double));
		if (tempInvA == NULL){
			printf("ERROR : not enough memory for array tempInvA in matrix \
			inversion\n");
			return EXIT_FAILURE;
		}
		/* temporary array in which every inverse values will be stored
		(even zeros) */

		/* Computing inverse of matrix A */

		int nnzInvA = 0; /* non zero elements of A^-1, this needs to be computed
		before generating the ctual arrays */
		for (int i = 0; i < problemSize; i++){
			// computing the columns of A^-1
			int *b = calloc(problemSize * sizeof(int));
			if (b == NULL){
				printf("ERROR : not enough memory for array b_i in matrix inversion\n");
				return EXIT_FAILURE;
			}
			b[i] = 1; // ith vector of the canonic basis
			double *ithColumnInvA = malloc(problemSize * sizeof(double)); /* ith
			column of A^-1, computed using UMF Pack */
			if (ithColumnInvA == NULL){
				printf("ERROR : not enough memory for array ithColumn in matrix\
				inversion\n");
			}
			if (umfSolve(problemSize, a, ja, ia, ithColumnInvA, b)){
				printf("ERROR : UMF Solve failed while inversing a matrix.\n");
				free(tempInvA); free(b); free(ithColumnInvA); // releasing memory
				return EXIT_FAILURE;
			}
			for (int j = 0; j < problemSize; j++){
				// iterating through the computed ith column of A^-1
				tempInvA[i * problemSize + j] = ithColumnInvA[j]; /* storing the jth
				element of the ith column into the temporary array */
				if (ithColumnInvA[j] != 0){
					nnzLA++; // increases the number of nnz elements
				}
			}
			free(b); free(ithColumnInvA); // freeing the memory
		}

		/* Formatting matrix A^-1 to CSR */

		// Memory allocations

		*invA = malloc(nnzInvA * sizeof(double)); /* allocation of memory of the
		array invA array of the matrix A^-1 */
		*invJa = malloc(nnzInvA * sizeof(int)); /* allocation of memory of the
		array invJa array of the matrix A^-1 */
		*invIa = malloc((problemSize + 1) * sizeof(int)); /* allocation of memory
		of the array invIa array of the matrix A^-1 */

		nnzInvA = 0; // back to zero, this is required to fill the arrays properly

		for (int i = 0; i < problemSize; i++){
			(*invIa)[i] = nnzInvA;
			for (int j = 0; j < problemSize; j++){
				/* iterating through all the elements of tempInvA and filling in the
				arrays */
				if (tempInvA[i + j * problemSize] != 0){
					(*invA)[nnz] = tempInvA[i + j * problemSize];
					(*invJa)[nnz] = j;
					nnzInvA++;
				}
			}
		}
		invIa[problemSize] = nnzInvA; //saving of nnzInvA
		free(tempInvA); // freeing memory

		return EXIT_SUCCESS;
	}
