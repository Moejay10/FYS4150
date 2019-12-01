#ifndef FUNCTION_H
#define FUNCTION_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace  std;
using namespace arma;



void forward_step(double, rowvec &, rowvec &, int);
void forward_euler(double, mat &, int, int);
void tridiag(double, rowvec &, int);
void backward_euler(double, mat &, int, int);
void crank_nicolson(double, mat &, int, int);
void g(mat &, int);
void analytic(mat &, int, int, vec, double);

#endif
