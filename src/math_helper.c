// Author : Denis Verstraeten
// Created on : 27/11/2017

/* this module is used for all calculations */

#include "math_helper.h"

int makeSuitableM(int *m){
	/* this funtions modifies m to a correct value of m so that the multi-grid
	method is usable. m needs to belong to 2^n + 1 so that the corresponding step
	is always divisible by 2. It returns the value of n.*/
	int isIncorrect = 1; // boolean value telling whether m is incorrect
	while(isIncorrect){
		for (int n = 1; n <= 20; n++){
			/* checks until 2^20 + 1 = 1048577, so practically much bigger than
			actual values of m */
			if(*m == pow(2, n) + 1){
				// if m is suitable
				return n;
			}
		}
		(*m)++; // increases m if not suitable
	}
	printf("ERROR : could not find a suitable value of m \n");
	return EXIT_FAILURE;
}

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

int fullSplitAMatrix(int m, int problemSize, double *a, int *ia, int *ja,
	double **la, int **ila, int **jla, double **ua, int **iua, int **jua,
	double **da, int **ida, int **jda){
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
	*jda = malloc(problemSize * sizeof(int)); /* allocation of memory for the jda
	array corresponding to the DA matrix */
	*ida = malloc((problemSize + 1) * sizeof(int)); /* allocation of memory of the
	ida array corresponding to the DA matrix */

	// checks whether the allocation was successful, returns an error if not

	if(*la == NULL){
		printf("ERROR : not enough memory to generate la array\n");
		return EXIT_FAILURE
	}

	if (*jla == NULL){
		printf("ERROR : not enough memory to generate jla array\n");
		free(*la);
		return EXIT_FAILURE;
	}

	if (*ila == NULL){
		printf("ERROR : not enough memory to generate ila array\n");
		free(*la); free(*jla);
		return EXIT_FAILURE;
	}

	if (*ua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila);
		return EXIT_FAILURE;
	}

	if (*jua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila); free(*ua);
		return EXIT_FAILURE;
	}

	if (*iua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila); free(*ua); free(*jua);
		return EXIT_FAILURE;
	}

	if (*da == NULL){
		printf("ERROR : not enough memory to generate da array\n");
		free(*la); free(*jla); free(*ila); free(*ua); free(*jua); free(*iua);
		return EXIT_FAILURE;
	}

	if (*jda == NULL){
		printf("ERROR : not enough memory to generate jda array\n");
		free(*la); free(*jla); free(*ila); free(*ua); free(*jua); free(*iua);
		free(*da);
		return EXIT_FAILURE;
	}

	if (*ida == NULL){
		printf("ERROR : not enough memory to generate ida array\n");
		free(*la); free(*jla); free(*ila); free(*ua); free(*jua); free(*iua);
		free(*da); free(*jda);
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
				(*jda)[i] = i;
				(*ida)[i] = i;
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
	(*ida)[problemSize] = problemSize; // saving of nnz of DA
	return EXIT_SUCCESS;
}

int notDiagonalSplitAMatrix(int m, int problemSize, double *a, int *ia, int *ja,
	double **la, int **ila, int **jla, double **ua, int **iua, int **jua){
	/* this funtions does the same thing as the previous one, except that it
	does not return the arrays corresponding to matrix D. Only the arrays linked
	to LA and UA are returned */

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

	// checks whether the allocation was successful, returns an error if not

	if(*la == NULL){
		printf("ERROR : not enough memory to generate la array\n");
		return EXIT_FAILURE
	}

	if (*jla == NULL){
		printf("ERROR : not enough memory to generate jla array\n");
		free(*la);
		return EXIT_FAILURE;
	}

	if (*ila == NULL){
		printf("ERROR : not enough memory to generate ila array\n");
		free(*la); free(*jla);
		return EXIT_FAILURE;
	}

	if (*ua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila);
		return EXIT_FAILURE;
	}

	if (*jua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila); free(*ua);
		return EXIT_FAILURE;
	}

	if (*iua == NULL){
		printf("ERROR : not enough memory to generate ua array\n");
		free(*la); free(*jla); free(*ila); free(*ua); free(*jua);
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
			printf("ERROR : not enough memory to generate x array\n");
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
			printf("ERROR : not enough memory to generate x array\n");
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

int S(int j){
	// this returns the sum 1+2+3+...+j
	return j * (j + 1) / 2;
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

		int maxNnzInvA = S(problemSize); /* this is the worst
		case scenario of filling in of the triangular matrix of dimension
		problemSize*problemSize */
		double *tempInvA = malloc(maxNnzInvA * sizeof(double));
		if (tempInvA == NULL){
			printf("ERROR : not enough memory for array tempInvA in matrix"
			" inversion\n");
			return EXIT_FAILURE;
		}
		/* temporary array in which every inverse values will be stored
		(even zeros) */

		/* Computing inverse of matrix A */

		int nnzInvA = 0; /* non zero elements of A^-1, this needs to be computed
		before generating the actual arrays */
		int nnzTempInvA = 0; // non zero elements of the temporary array
		for (int i = 0; i < problemSize; i++){
			// computing the columns of A^-1
			double *b = malloc(problemSize * sizeof(double));
			if (b == NULL){
				printf("ERROR : not enough memory for array b_i in matrix inversion\n");
				free(tempInvA);
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
					printf("ERROR : upper triangular matrix solving failed to inverse"
					" matrix\n");
					free(b); free(tempInvA);
					return EXIT_FAILURE;
				}
				for (int j = 0; j <= i; j++){
					// iterating through the computed ith column of A^-1
					tempInvA[nnzTempInvA] = ithColumnInvA[j]; /* storing the jth
					element of the ith column into the temporary array */
					nnzTempInvA++; /* increases the number of non zero element of
					temporary array */
					if (ithColumnInvA[j] != 0){
						nnzInvA++; // increases the number of nnz elements
					}
				}
			} else {
				// the matrix to inverse is a lower triangular one
				if(lowerTriangularSolver(problemSize, a, ja, ia, &ithColumnInvA, b)){
					printf("ERROR : lower triangular matrix solving failed to inverse"
					" matrix\n");
					free(b); free(tempInvA);:;lk
					return EXIT_FAILURE;
				}
				for (int j = i; j < problemSize; j++){
					// iterating through the computed ith column of A^-1
					tempInvA[nnzTempInvA] = ithColumnInvA[j]; /* storing the jth
					element of the ith column into the temporary array */
					if (ithColumnInvA[j] != 0){
						nnzInvA++; // increases the number of nnz elements
						nnzTempInvA++; /* increases the number of non zero element of
						temporary array */
					}
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

		// Memory check

		if (*invA == NULL){
			printf("ERROR : not enough memory to generate array invA\n");
			return EXIT_FAILURE;
		}

		if (*invJa == NULL){
			printf("ERROR : not enough memory to generate array invJa\n");
			free(*invA);
			return EXIT_FAILURE;
		}

		if (*invIa == NULL){
			printf("ERROR : not enough memory to generate array invIa\n");
			free(*invA); free(*invJa);
			return EXIT_FAILURE;
		}

		nnzInvA = 0; // back to zero, this is required to fill the arrays properly

		if (UP){
			for (int i = 0; i < problemSize; i++){
				(*invIa)[i] = nnzInvA;
				for (int j = i; j < problemSize; j++){
					/* iterating through all the elements of tempInvA and filling in the
					arrays */
					if (tempInvA[i + S(j)] != 0){
						(*invA)[nnzInvA] = tempInvA[i + S(j)];
						(*invJa)[nnzInvA] = j;
						nnzInvA++;
					}
				}
			}
		} else {
			for (int i = 0; i < problemSize; i++){
				(*invIa)[i] = nnzInvA;
				for (int j = 0; j <= i; j++){
					/* iterating through all the elements of tempInvA and filling in the
					arrays */
					int tempIndex = i + S(problemSize - 1) - S(problemSize - 1 - j);
					/* index of the tempInvA corresponding to i and j */
					if (tempInvA[tempIndex] != 0){
						(*invA)[nnzInvA] = tempInvA[tempIndex];
						(*invJa)[nnzInvA] = j;
						nnzInvA++;
					}
				}
			}
		}
		(*invIa)[problemSize] = nnzInvA; //saving of nnzInvA
		free(tempInvA); // freeing memory

		return EXIT_SUCCESS;
	}

int matrixVectorMultCSR(int m, double *a, int *ja, int *ia,
	double *vector, double **newVector){
		/* this function computes the product A*vector with A a matrix and
		vector a vector of matching dimensions
		dim(newVector) = mx1
		dim(A) = mxn
		dim(vector) = nx1 */

		for (int i = 0; i < m; i++){
			// iterating through the lines of A and vector to compute the product
			double sum = 0; // value of the ith component of the vector
			for (int j = ia[i]; j < ia[i + 1]; j++){
				// iterating through the non zero elements of the ith line of A
				sum += a[j] * vector[ja[j]];
			}
			(*newVector)[i] = sum;
			sum = 0;
		}
		return EXIT_SUCCESS;
	}

double normVector(int problemSize, double *vector){
	// this function returns the Euclidean norm of a vector
	double squaredNorm = 0; // sum of the squared elements of the vector
	for (int i = 0; i < problemSize; i++){
		// iterating through the elements of vector
		squaredNorm += square(vector[i]);
	}
	return sqrt(squaredNorm); // returns the square root of the squared norm
}

int vectorDifference(int problemSize, double *vector1, double *vector2,
	double **difference){
	/* this functions returns the difference between the two vectors vector 1 and
	vector2 */

	for (int i = 0; i < problemSize; i++){
		// iterating through the elements of vector1, vector 2 and difference
		(*difference)[i] = vector1[i] - vector2[i];
	}

	return EXIT_SUCCESS;
}

int vectorAddition(int problemSize, double *vector1, double *vector2,
	double **addition){
		// this function returns the sum of the two vectors vector1 and vector2

		for (int i = 0; i < problemSize; i++){
			// iterating through the elements of vector1, vector 2 and difference
			(*addition)[i] = vector1[i] + vector2[i];
		}

		return EXIT_SUCCESS;
	}

int getResidue(int problemSize, double *b, double *a, int *ja, int *ia,
	double *u, double **residue){
		/* this function returns the residue of the current approximation by
		computing r_m = b - A*u_m */

		// Memory allocation

		double *product = malloc(problemSize * sizeof(double)); /* array containing
		the result of the matrix-vector product */

		// Memory check

		if (product == NULL){
			printf("ERROR : not enough memory to generate product array\n");
			return EXIT_FAILURE;
		}

		// Computing A*u_m

		if (matrixVectorMultCSR(problemSize, a, ja, ia, u, &product)){
			printf("ERROR : could not compute matrix-vector product\n");
			return EXIT_FAILURE;
		}

		// Computing b - A*u_m

		if (vectorDifference(problemSize, b, product, residue)){
			printf("ERROR : could not compute vector difference\n");
			free(product);
			return EXIT_FAILURE;
		}

		free(product); // freeing memory

		return EXIT_SUCCESS;

	}
