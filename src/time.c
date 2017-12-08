// Author : Denis Verstraeten
// Created on : 08/12/2017

/* this module is concerned with functions about the time */

#include <time.h>

double timer(void){
  // this functions returns the time since the beginning of the execution
  return (double) clock() / CLOCKS_PER_SEC;
}
