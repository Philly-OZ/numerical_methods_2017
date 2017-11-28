// Author : Denis Verstraeten
// Created on : 28/11/2017

/* this module will be used to plot the results using the GNUPLOT library */

#include "visual.h"

int plot(int m, double step, double *T, double *dirichletCond){
  // this function will plot the vector of solution in a 2D heat map

  /* Connecting to the files and saving the data */

  FILE* gnuplot = popen("gnuplot -persist", "w"); // opens the GNUPLOT file
  FILE* data = fopen("data.txt", "w"); /* opens the file in which the data
  will be temporarily stored before being displayed */

  if (gnuplot == NULL || data == NULL){
    printf("ERROR : cannot plot result\n");
    return EXIT_FAILURE;
  }

  double x, y; // cartesian coordinates
  int equationNumber; // number of the equation

  for (int i = 0; i < m; i++){
    y = i * step;
    for (int j = 0; j < m - 1; j++){
      // iterating through the solution and storing them into a chart
      x = j * step;
      equationNumber = i * (m - 1) + j;
      fprintf(data, "%f %f %f\n", x, y, T[equationNumber]); // saving the data
    }
    x = (m - 1) * step; // points on the Eastt edge (Dirichlet condition)
    fprintf(data, "%f %f %f\n", x, y, dirichletCond[equationNumber]);
  }

  fclose(data); // terminates the saving of data

  /* Formatting the chart */

  fprintf(gnuplot, "set term png font '/Library/Fonts/Arial.ttf' 14\n");
  fprintf(gnuplot, "set size square\n");
  fprintf(gnuplot, "set output 'heat_map.png'\n");
  fprintf(gnuplot, "set view map\n");
  fprintf(gnuplot, "set dgrid3d\n");
  fprintf(gnuplot, "set title 'Stationary solution fo the heat equation'\n");
  fprintf(gnuplot, "set xlabel 'x [m]'\n");
  fprintf(gnuplot, "set xrange [0:0.2]\n");
  fprintf(gnuplot, "set ylabel 'y [m]'\n");
  fprintf(gnuplot, "set yrange [0:0.2]\n");
  //fprintf(gnuplot, "set cbrange [0:5]\n");
  fprintf(gnuplot, "set pm3d interpolate 10,10\n");
  fprintf(gnuplot, "splot 'data.txt' using 1:2:3 with pm3d\n");

  pclose(gnuplot); // terminates the plot

  remove("data.txt"); // removes the temporary saved data
  system("mv heat_map.png graphics"); // moves the png to the output directory

  return EXIT_SUCCESS;
}
