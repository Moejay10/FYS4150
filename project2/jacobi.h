#ifndef JACOBI_H
#define	JACOBI_H

#include <iostream>
#include <cmath>
#include <armadillo>
#include <tuple>

using namespace std;
using namespace arma;

tuple <double,int,int> Find_max_nondiagonal(mat);
mat Jacobi(mat, double, int);

#endif /* JACOBI_H */
