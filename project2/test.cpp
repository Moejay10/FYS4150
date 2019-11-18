#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "jacobi.h"

TEST_CASE( "Testing the maximal non-diagonal finder"){
  
  // Create a test matrix
  int n = 5;
  mat A = zeros<mat>(n,n);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      A(i,j) = i + 2*j;
    }
  }

  // Call the Find_max_nondiagonal-function
  double max = get<0>(Find_max_nondiagonal(A));
  int k = get<1>(Find_max_nondiagonal(A));
  int l = get<2>(Find_max_nondiagonal(A));

  // Check the values
  CHECK( max == 11.0 );
  CHECK( k == 3 );
  CHECK( l == 4 );
}

TEST_CASE( "Testing orthogonality of the eigenvectors"){

  // Create the test matrix
  int n = 3;
  vec a  = linspace<vec>(1,n,n);
  mat A = toeplitz(a);
  
  // Declare vectors and matrices for the eigenvalues and eigenvectors before and after the transformation
  vec eigval_before; vec eigval_after; mat eigvec_before; mat eigvec_after;

  // Use armadillo to obtain the eigenpairs (before)
  eig_sym(eigval_before, eigvec_before, A);

  // Perform the transformation
  mat  B = Jacobi(A, 1e-8, 1000);

  // Use armadillo to obtain the eigenpairs (after)
  eig_sym(eigval_after, eigvec_after, A);

  // Calculate the dot product of the different eigenvectors with each other
  double dot1_before = eigvec_before(0,0)*eigvec_before(0,1) + eigvec_before(1,0)*eigvec_before(1,1) +  eigvec_before(2,0)*eigvec_before(2,1); // v1 * v2 
  double dot2_before = eigvec_before(0,0)*eigvec_before(0,2) + eigvec_before(1,0)*eigvec_before(1,2) +  eigvec_before(2,0)*eigvec_before(2,2); // v1 * v3
  double dot3_before = eigvec_before(0,1)*eigvec_before(0,2) + eigvec_before(1,1)*eigvec_before(1,2) +  eigvec_before(2,1)*eigvec_before(2,2); // v2 * v3

  double dot1_after = eigvec_after(0,0)*eigvec_after(0,1) + eigvec_after(1,0)*eigvec_after(1,1) +  eigvec_after(2,0)*eigvec_after(2,1);        // v1 * v2
  double dot2_after = eigvec_after(0,0)*eigvec_after(0,2) + eigvec_after(1,0)*eigvec_after(1,2) +  eigvec_after(2,0)*eigvec_after(2,2);        // v1 * v3
  double dot3_after = eigvec_after(0,1)*eigvec_after(0,2) + eigvec_after(1,1)*eigvec_after(1,2) +  eigvec_after(2,1)*eigvec_after(2,2);        // v2 * v3

  CHECK( dot1_before == Approx(0).margin(1e-10) );
  CHECK( dot2_before == Approx(0).margin(1e-10) );
  CHECK( dot3_before == Approx(0).margin(1e-10) );
  CHECK( dot1_after == Approx(0).margin(1e-10) );
  CHECK( dot2_after == Approx(0).margin(1e-10) );
  CHECK( dot3_after == Approx(0).margin(1e-10) );
}
