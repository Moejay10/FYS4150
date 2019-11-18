#include "jacobi.h"

using namespace chrono;

int main(int argc, char *argv[]){
  // Declare the number of steps and the maximum number of iterations
  int N; int max_iterations;
    // Read the number of steps and the number of iterations
    if( argc <= 2 ){
      cout << "Error: " << argv[0] << " reads the number of steps and the maximal number of iterations" << endl;
      exit(1);
    }
    else{
      N = atoi(argv[1]);
      max_iterations = atoi(argv[2]);
    }

  // Define the step size
  double h = 1.0/N;
  double hh = h*h;
  
  // Define the matrix
  mat A = zeros<mat>(N-1,N-1);
  // Set up the matrix elements
  double d = 2.0/hh; double a = -1.0/hh; 
  A(N-2,N-2) = d;
  for (int i = 0; i <= N-3; i++){
    A(i,i+1) = a;
    A(i,i) = d;
    A(i+1,i) = a;
  }

  // Define a vector for the analytical eigenvalues
  vec eigenvalues_analytic(N-1);
  // Define pi
  double pi = 4*atan(1);
  // Calculate the analytical eigenvalues
  for (int i = 0; i <= N-2; i++){ eigenvalues_analytic(i) = d + 2*a*cos((i+1)*pi/N);}

  // Finding eigenvalues using Armadillo
  vec eigenvalues_arma(N-1);
  // Time the algorithm
  auto start = high_resolution_clock::now();
  eig_sym(eigenvalues_arma, A);
  auto finish = high_resolution_clock::now();
  duration<double> time_used = finish - start;

  // Print the time
  printf("Time used by Armadillo: %.3e\n",time_used.count());

  // Diagonalise the matrix using Jacobi
  A = Jacobi(A, 1e-8, max_iterations);
  
  // Define a vector for the numerical eigenvalues
  vec eigenvalues(N-1);
  // Collect the eigenvalues from the diagonal matrix
  for (int i = 0; i <= N-2; i++){
    eigenvalues(i) = A(i,i);
  }
  // Sort the eigenvalues
  eigenvalues = sort(eigenvalues);

  // Print the numerical and analytical eigenvalues and the relative error
  printf(" Analytical   Numerical   Relative error\n");
  for (int i = 0; i <= N-2; i++){
    double error = fabs((eigenvalues_analytic(i) - eigenvalues(i))/eigenvalues_analytic(i));
    printf("%10.5f   %10.5f   %.3e\n", eigenvalues_analytic(i), eigenvalues(i), error);
  }

  return 0;
}
