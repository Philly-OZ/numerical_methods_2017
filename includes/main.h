// Libraries

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// Function declarations

int generate_problem(int, double, double, int *, int **, int **, double **,
double **, double **);
void generateOptions(int *, int *, int *, int *, int *, int *, int *,
  int, char **);
void printAArrays(double *, int *, int *, int);
void printLinearSystemArrays(double *, int *, int *, double *, int);
int umfSolve(int, double *, int *j, int *, double *, double *);
int plot(int, double, double *, double *);
double timer(void);
time_t time(int);
int sgsSolve(double *, int *, int *, double **, double *, double, int, int,
   int, int, int);
