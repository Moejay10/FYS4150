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
#include <omp.h>

using namespace  std;
using namespace arma;


// One Dimensional case functions
void forward_step(double, rowvec &, rowvec &, int, bool);
void forward_step_2dim(double, mat &, mat &, int);
void forward_euler(double, mat &, int, int);
void forward_euler_2dim(double, cube &, int, int);
void tridiagSolver(rowvec &, rowvec, double, int, bool);
void backward_euler(double, mat &, int, int);
void crank_nicolson(double, mat &, int, int);
void analytic(mat &, int, int, int);

// Two Dimensional case functions
void analytic_2D(mat &, int, double);
int JacobiSolver(mat &, double, double, double, int);

// Lithosphere functions
void Lithosphere(int, double, double, double, int);
void Heat(mat &, int, int);
void Qzones(vec &, double, double, int);
void Decay(vec &, int T, double Q_s);


#endif
