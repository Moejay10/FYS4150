#include <iostream>
#include <cmath>
#include <armadillo>
#include <tuple>
#include <chrono>

#include "jacobi.h"

using namespace chrono;

// Find the largest non-diagonal element

tuple <double,int,int> Find_max_nondiagonal(mat A){
  // Declare the maximal non-diagonal element, and its indices k and l
  double max; int k, l;
  // Define the dimension of the matrix
  int N = A.n_cols;
  // Find max, k and l
  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      double a_ij = fabs(A(i,j));
      if (a_ij > max){
	max = a_ij; k = i; l = j;
      }
    }
  }
return make_tuple(max, k, l);
}


// Perform a Jacobi rotation, diagonalising A and finding the eigenvalues

mat Jacobi(mat A, double epsilon, int max_iterations){
  // Collect the maximal non-diagonal element and its indices from Find_max_nondiagonal()
  double max_nondiagonal = get<0>(Find_max_nondiagonal(A));
  int k = get<1>(Find_max_nondiagonal(A));
  int l = get<2>(Find_max_nondiagonal(A)); 
  
  // Time the algorithm
  auto start = high_resolution_clock::now();

  // Perform the Jacobi rotation until the maximal non-diagonal element reaches the desired tolerance, or after a maximum number of iterations
  int iterations = 0;
    while (max_nondiagonal > epsilon && iterations < max_iterations){

    // Declare tau, tan, sin and cos
    double tau,t,s,c;

    // Prevent divison by zero
    if (A(k,l) == 0){
      s = 0.0; c = 1.0;
    }
    else{
      // Calculate tau
      tau = (A(l,l) - A(k,k))/(2*A(k,l));
      // Find the smallest solution for t
      if (tau > 0){
	t = 1.0/(tau + sqrt(1.0 + tau*tau));
      }
      else{
	t = -1.0/(- tau + sqrt(1.0 + tau*tau));
      }
      // Calculate c and s
      c = 1.0/sqrt(1.0 + t*t); 
      s = c*t;
    }

    // Perform the similarity transformation
    double a_kk = A(k,k); double a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = c*c*a_ll + 2*c*s*A(k,l) + s*s*a_kk;
    A(k,l) = 0.0; A(l,k) = 0.0;
    for (int i = 0; i < A.n_cols; i++){
      if (i != k && i != l){
	double a_ik = A(i,k); double a_il = A(i,l);
	A(i,k) = c*a_ik - s*a_il;
	A(i,l) = c*a_il + s*a_ik;
	A(k,i) = A(i,k);
	A(l,i) = A(i,l);
      }
}

    // Update the maximal non-diagonal element and its indices
    max_nondiagonal = get<0>(Find_max_nondiagonal(A));
    k = get<1>(Find_max_nondiagonal(A));
    l = get<2>(Find_max_nondiagonal(A));

    iterations++;
    }

  auto finish = high_resolution_clock::now();
  duration<double> time_used = finish - start;
  
  // Print the number of iterations and time used
  printf("Number of iterations: %d\n", iterations);
  printf("Time used by Jacobi: %.3e\n", time_used.count());
  return A;
}
