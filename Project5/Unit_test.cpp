#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "Functions.h"

TEST_CASE( "Checking the Numerical vs Analytical results in 1 dimension." ){

// Define parameters for the test case in one dimension
  double PI = 4*atan(1);

  // Step length in position x
  double dx = 0.1;

  // Step length in time
  double dt = 0.5*dx*dx;


  // Defining alpha
  double alpha = dt/dx/dx;

  // Number of integration points along x-axis (inner points only)
  int N = int(1.0/(dx));

  // Number of time steps till final time
  int T = int(1.0/dt);

  // Defining u
  mat u = zeros<mat>(T,N+2);
  mat u_analytic1 = zeros<mat>(T,N+2);


  // Implement boundaries rigidly
  // Boundary condition for first endpoint is already zero, so just need to
  // implement boundary for last end point
  for (int t = 0; t < T; t++){
    u(t,N+1) = 1.0;
  }



  // Analytical results for diffusion equation in one dimension
  int Max = 1000; // what is defined as numerical "infinity".
  analytic(u_analytic1, N, T, Max);

  // Numerical methods in one dimension
  mat uf = u;
  mat ub = u;
  mat uc = u;
  forward_euler(alpha, uf, N, T);
  backward_euler(alpha,ub,N,T);
  crank_nicolson(alpha, uc, N, T);

  // Check the values in one dimension
  double tol = 1e-1;
  for (int t = 0; t < T; t++){
    for (int i = 0; i < N+1; i++){
  CHECK( fabs(uf(t,i) - u_analytic1(t,i)) < 2*tol);
  CHECK( fabs(ub(t,i) - u_analytic1(t,i)) < tol);
  CHECK( fabs(uc(t,i) - u_analytic1(t,i)) < tol);

    }
  }

}

TEST_CASE( "Checking the Numerical vs Analytical results in 2 dimensions for t = 1.0." ){


// Define parameters for the test case in two dimensions
  double PI = 4*atan(1);

  // Step length in position x
  double dx = 0.1;

// Step length in time
  double dt = 0.2*dx*dx;

  // Defining alpha
  double alpha = dt/dx/dx;

  // Number of integration points along x-axis (inner points only)
  int N = int(1.0/(dx));

// Number of time steps till final time
  int T = int(1/dt);

  // Analytical solution to the diffusion equation in two dimensions
  // Defining u
  mat u_analytic2 = zeros<mat>(N+2,N+2);
  double temp;
  for (int t = 0; t < T; t++){
    temp = t*dt;
    analytic_2D(u_analytic2, N, temp);
  }

  // Explicit method in two dimensions
  cube uf2 = zeros<cube>(N+2, N+2, T);
  //Boundary conditions
  for (int t = 0; t < T; t++){
    for (int n = 0; n < N+2; n++){
      uf2(n,N+1,t) = 0.0;
      uf2(N+1,n,t) = 0.0;
    }
  }
  // Initial conditions
    for (double i = 1; i < N+1; i++) {
      for (double j =1; j < N+1; j++) {
        uf2(i,j,0) = sin(i*PI/(double)N)*sin(j*PI/(double)N);
      }
    }
    forward_euler_2dim(alpha,uf2,N,T);

  // Implicit method in two dimensions
    mat u_implicit = zeros<mat>(N+2, N+2);


    // Implement boundaries rigidly
    // Boundary condition for endpoints is already zero
    //Boundary conditions
    for (int n = 0; n < N+2; n++){
      u_implicit(n, N+1) = 0.0;
      u_implicit(N+1, n) = 0.0;
    }
    // Initial conditions
      for (double i = 1; i < N+1; i++) {
        for (double j =1; j < N+1; j++) {
          u_implicit(i,j) = sin(i*PI/(double)N)*sin(j*PI/(double)N);
        }
      }

    double tolerance = 1e-10; // tolerance for convergence in Jacobi method
    int maxiter = 10000; // Max no. of iterations each time step in Jacobi Method

    int itcount = JacobiSolver(u_implicit, dx, dt, tolerance, maxiter);


    // Check the values in two dimensions for t = 1.0
    double tol = 1e-6;
    for (int i = 0; i < N+1; i++){
      for (int j = 0; j < N+1; j++){
    CHECK( fabs(u_implicit(i,j) - u_analytic2(i,j)) < tol);
    CHECK( fabs(uf2(i,j,T-1) - u_analytic2(i,j)) < tol);
    }
  }

}
