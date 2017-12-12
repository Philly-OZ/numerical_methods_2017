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

int lowerTriangularSolver(int problemSize, double *a, int *ja, int *ia,
	double **x, double *b){
		/* this function solves the system Ax=b with A a lower triangular matrix
		stored in CSR format */
		*x = malloc(problemSize * sizeof(double)); /* memory allocation for solution
		array */
		if (*x == NULL){
			printf("ERROR : not enough memory to inverse the lower triangular \
			matrix.\n");
			return EXIT_FAILURE;
		}
		for (int i = 0; i < problemSize; i++){
			// iterating through the lines of the matrix A
			double sum = 0; // sum_{j<i} a[i,j]*x[j]
			for (int j = ia[i]; j < ia[i + 1] - 1; j++){
				// iterating through the non zeros elements of the ith line of matrix A
				sum += a[j] * (*x)[ja[j]];
			}
			(*x)[i] = (1.0 / a[ia[i + 1] - 1]) * (b[i] - sum);
			if ((*x)[i] < 1e-6 && (*x)[i] > - 1e-6){
				// to avoid a 0 being saved as 0.0000000.....
				(*x)[i] = 0;
			}
			sum = 0;
		}
		return EXIT_SUCCESS;
	}

int upperTriangularSolver(int problemSize, double *a, int *ja, int *ia,
	double **x, double *b){
		/* this function solves the system Ax=b with A a upper triangular matrix
		stored in CSR format */
		*x = malloc(problemSize * sizeof(double)); /* memory allocation for solution
		array */
		if (*x == NULL){
			printf("ERROR : not enough memory to inverse the upper triangular \
			matrix.\n");
			return EXIT_FAILURE;
		}
		for (int i = problemSize - 1; i >= 0; i--){
			// iterating through the lines of the matrix A
			double sum = 0; // sum_{j<i} a[i,j]*x[j]
			for (int j = ia[i] + 1; j < ia[i + 1]; j++){
				// iterating through the non zeros elements of the ith line of matrix A
				sum += a[j] * (*x)[ja[j]];
			}
			(*x)[i] = (1.0 / a[ia[i]]) * (b[i] - sum);
			if ((*x)[i] < 1e-6 && (*x)[i] > - 1e-6){
				// to avoid a 0 being saved as 0.0000000.....
				(*x)[i] = 0;
			}
			sum = 0;
		}
		return EXIT_SUCCESS;
	}

int inverseMatrix(int problemSize, double **invA, int **invJa, int **invIa,
	double *a, int *ja, int *ia, int UP){
		/* this functions returns the inverse of the matrix A stored in CSR in
		a, ja and ia as CSR arrays invA, invJa, invIa. It uses the UMF Pack solver
		since this will actually only deals with triangular matrix.
		It solves the system Ax=b <=> x=A\b with b a vector of the canonic basis
		=> solving the system b_i gives the ith column of A^-1.
		UP parameter is a boolean value telling whether the matrix is upper or lower
		triangular */

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
			double *b = malloc(problemSize * sizeof(double));
			if (b == NULL){
				printf("ERROR : not enough memory for array b_i in matrix inversion\n");
				return EXIT_FAILURE;
			}
			for (int j = 0; j < problemSize; j++){
				// forming the ith vector of the canonic basis
				if (j != i){
					b[j] = 0.0;
				} else {
					b[j] = 1.0;
				}
			}

			double *ithColumnInvA; /* ith column of A^-1, computed using triangular
			solvers */

			if (UP){
				// the matrix to inverse is a upper triangular one
				if (upperTriangularSolver(problemSize, a, ja, ia, &ithColumnInvA, b)){
					printf("ERROR : upper triangular matrix solving failed to inverse\
					matrix\n");
					free(b);
					return EXIT_FAILURE;
				}
			} else {
				// the matrix to inverse is a lower triangular one
				if(lowerTriangularSolver(problemSize, a, ja, ia, &ithColumnInvA, b)){
					printf("ERROR : lower triangular matrix solving failed to inverse\
					matrix\n");
					free(b);
					return EXIT_FAILURE;
				}
			}

			for (int j = 0; j < problemSize; j++){
				// iterating through the computed ith column of A^-1
				tempInvA[i * problemSize + j] = ithColumnInvA[j]; /* storing the jth
				element of the ith column into the temporary array */
				if (ithColumnInvA[j] != 0){
					nnzInvA++; // increases the number of nnz elements
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
					(*invA)[nnzInvA] = tempInvA[i + j * problemSize];
					(*invJa)[nnzInvA] = j;
					nnzInvA++;
				}
			}
		}
		(*invIa)[problemSize] = nnzInvA; //saving of nnzInvA
		free(tempInvA); // freeing memory

		return EXIT_SUCCESS;
	}

double getIJElementCSR(int i, int j, double *a, int *ja, int *ia){
	// this funtions returns the value A[i,j] of a matrix A stored in CSR format
	for(int k = ia[i]; k < ia[i + 1]; k++){
		// iterating through the non zero elements of the ith line of A
		if (j == ja[k]){
			return a[k];
		}
	}
	return 0.0;
}

int makePreconditionner(int problemSize, double *da, double *invLa,
	double *invUa, int *invJla, int *invJua, int *invIla, int *invIua,
	double **prec, int **jPrec, int **iPrec){
		/* this funtion computes the product of the matrixes LA^-1, UA^-1 and DA to
		obtain the preconditionning B^-1 = UA^-1 * DA * LA^-1. It stores the result
		in CSR format */

		double *tempPrec = malloc(square(problemSize) * sizeof(double));
		/* temporary array in which every values of B^-1 will be stored
		(even zeros) */
		if(tempPrec == NULL){
			printf("ERROR : not enough memory to generate preconditionning matrix\n");
			return EXIT_FAILURE;
		}

		// Computing the product of the matrixes

		int nnzPrec = 0; // number of non zero elements of the preconditionner
		for (int i = 0; i < problemSize; i++){
			// iterating through the lines of B^-1
			for (int j = 0; j < problemSize; j++){
				// iterating through the elements of ith line of B^-1
				double sum = 0; // value that will be saved as B^-1[i,j]
				for (int k = 0; k < problemSize; k++){
					/* iterating through the ith line of UA^-1, the jth column of LA^-1
					and DA to compute the product */
					sum += getIJElementCSR(i, k, invUa, invJua, invIua) *
						getIJElementCSR(k, j, invLa, invJla, invIla) * da[k];
				}
				if (sum != 0){
					nnzPrec ++;
				}
				tempPrec[i * problemSize + j] = sum; /* saving the computed value in the
				temporary array */
			}
		}

		/* Formatting the result to CSR */

		// Memory allocations

		*prec = malloc(nnzPrec * sizeof(double));
		*jPrec = malloc(nnzPrec * sizeof(double));
		*iPrec = malloc((problemSize + 1) * sizeof(double)); /* allocation of memory
		for the array in which B^-1 will be stored */

		if (*prec == NULL || *jPrec == NULL || *iPrec == NULL){
			printf("ERROR :  not enough memory to store preconditionner\n");
			free(tempPrec);
			return EXIT_FAILURE;
		}

		// Filling the arrays

		nnzPrec = 0; // back to 0, it will be used to fill in the arrays

		for (int i = 0; i < problemSize; i++){
			//iterating through the lines of tempPrec
			(*iPrec)[i] = nnzPrec;
			for (int j = 0; j < problemSize; j++){
				// iterating through the elements of the ith line of tempPrec
				if (tempPrec[i * problemSize + j] != 0.0){
					// non zero element, needs to be saved in the arrays
					(*prec)[nnzPrec] = tempPrec[i * problemSize + j];
					(*jPrec)[nnzPrec] = j;
					nnzPrec++;
				}
			}
		}

		(*iPrec)[problemSize] = nnzPrec; //saving of the number of non zero elements

		return EXIT_SUCCESS;
	}
