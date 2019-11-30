#include "Functions.h"
#include "omp.h"


void forward_euler(double dx, mat &u, vec xlist, vec tlist){
  double C = dt/(dx*dx)
  double dt = 0.5*dx*dx //Stability criterion
  vec u1 = u(1) //Initial condition, zeros
  for(int j = 1; j < tlist.n_elem;  j++){
    for (int i = 1; i < xlist.n_elem-1; i++){
      u(j,i) = C*(u1(i+1) - 2*u1(i) + u1(i-1)) + u1(i)
    }
    // Boundaries need not be updated, they are already zeros by def. of u-matrix
    u1 = u(j)
  }
}

void backward_euler(double dx, mat &u, vec xlist, vec tlist){

}
void backward_euler(int n, int tsteps, double delta_x, double alpha)
{
   double a, b, c;
   vec u(n+1); // This is u  of Au = y
   vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step

   // Initial conditions
   for (int i = 1; i < n; i++) {
      y(i) = u(i) = func(delta_x*i);
   }
   // Boundary conditions (zero here)
   y(n) = u(n) = u(0) = y(0);
   // Matrix A, only constants
   a = c = - alpha;
   b = 1 + 2*alpha;
   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      tridag(a, b, c, y, u, n+1);
      // boundary conditions
      u(0) = 0;
      u(n) = 0;
      // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
	 y(i) = u(i);
      }
      //  You may consider printing the solution at regular time intervals
      ....   // print statements
   }  // end time iteration
   ...
}
