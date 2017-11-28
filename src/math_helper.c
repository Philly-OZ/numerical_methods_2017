// Author : Denis Verstraeten
// Created on : 27/11/2017

/* this module will be used for all calculations */

#include "math_helper.h"

int factorisation(int problemSize, double *a, int *ja, int *ia, void **Numeric){
  // this factorises a matrix to solve it
	int statut;
	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
	void *Symbolic;
	statut = umfpack_di_symbolic(problemSize, problemSize, ia, ja, a, &Symbolic,\
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

int resolution(int *ia, int *ja, double *a, double *x, double *b, \
  void **Numeric){
  // this will solve the linear system Ax=b with the matrix A already factorised
	int statut;
	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
	statut = umfpack_di_solve(UMFPACK_At, ia, ja, a, x, b, *Numeric, Control,\
     Info);
	if (statut < 0)
	{
		printf("ERROR : UMF_PACK_SOLVE FAILED\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

int umfSolve(int problemSize, double *a, int *ja, int *ia, double *x, \
  double *b){
    // this is called to actually solve the linear problem Ax=b
    void *Numeric = NULL; // variable required for the factorisation
    if (factorisation(problemSize, a, ja, ia, &Numeric) || \
  resolution(ia, ja, a, x, b, &Numeric)){
      return EXIT_FAILURE;
    } else {
      return EXIT_SUCCESS;
    }
  }
