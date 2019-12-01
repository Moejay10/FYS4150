#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include "omp.h"
#include "Functions.h"
using namespace  std;
using namespace arma;



int main(int argc, char* argv[]){

  cout << "\n" << "Project 5: Partial Differential Equations " << endl;
  cout << "\n" << "Which dimension of the diffusion equation do you want to run?: " << endl;
  cout << "\n" << "Project 5 - 1-dim: " <<  "Write 1 " << endl;
  cout << "\n" << "Project 5 - 2-dim: " <<  "Write 2 " << endl;



  cout << "\n" << "Write here " << endl;
  string Dimension;
  cin >> Dimension;

  // 1-Dimensional Diffusion Equation - Analytical vs Numerical Results for a 2x2 Lattice

  if (Dimension == "1"){



  cout << "\n" << "Choose step size for Delta x: " << endl;
  cout << "\n" << "Delta x = 0.1: " <<  "Write 0.1 " << endl;
  cout << "\n" << "Delta x= 0.01 " <<  "Write 0.01 " << endl;

  // Step length in position x
  double dx;
  cin >> dx;

  // Step length in time
  double dt = 0.5*dx*dx;


  // Defining alpha
  double alpha = dt/dx/dx;

  // Number of integration points along x-axis (inner points only)
  int N = int(1.0/(dx));

  // Number of time steps till final time
  int T = N;

  // Defining u
  mat u = zeros<mat>(T,N+2);
  mat u_analytic = zeros<mat>(T,N+2);



  vec x = linspace<vec>(0,1,N+2);


  // Implement boundaries rigidly
  // Boundary condition for first endpoint is already zero, so just need to
  // implement boundary for last end point
  for (int t = 0; t < T; t++){
    u(t,N+1) = 1.0;
    u_analytic(t,N+1) = 1.0;
  }

  // Initial codition
  g(u, N);
  g(u_analytic, N);

  analytic(u_analytic, N, T, x, dt);

//  u = u.t();
//  u_analytic = u_analytic.t();

  cout << "\n" << "Which method do you want to use to solve the diffusion equation?: " << endl;
  cout << "\n" << "Forward Euler- Explicit Scheme: " <<  "Write FE " << endl;
  cout << "\n" << "Backward Euler- Implicit Scheme: " <<  "Write BE " << endl;
  cout << "\n" << "Crank-Nicolson- Implicit Scheme: " <<  "Write CN " << endl;


  cout << "\n" << "Write here " << endl;
  string Method;
  cin >> Method;

if (Method == "FE"){
  forward_euler(alpha, u, N, T);
  cout << u << endl;

}
if (Method == "BE"){
  backward_euler(alpha,u,N,T);
  cout << u << endl;
}
if (Method == "CN"){
  crank_nicolson(alpha, u, N, T);
  cout << u << endl;
}





  }


  if (Dimension == "2"){


    cout << "\n" << "In progress: " << endl;
  }


  return 0;

}
