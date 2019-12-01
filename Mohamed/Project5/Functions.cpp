#include "Functions.h"
#include "omp.h"


// Code for solving the 1+1 dimensional diffusion equation
// du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
// with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0

double pi = 4*atan(1); // This is the constant pi

void forward_step(double alpha, rowvec &u, rowvec &uPrev, int N){
    /*
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    */

    for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (1.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
    //return u;
}

void forward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements the forward Euler sheme, results saved to
    array u
    */

    // Skip boundary elements
    rowvec u_temp_row;
    rowvec u_temp_row1;
    for (int t = 1; t < T; t++){

        u_temp_row = u.row(t);
        u_temp_row1 = u.row(t-1);

        forward_step(alpha,u_temp_row,u_temp_row1,N);
        u.row(t) = u_temp_row;

      }

}

void tridiag(double alpha, rowvec &u, int N){
    /*
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    */

    vec d(N);
    d.fill(1+2*alpha);
    vec b(N-1);
    b.fill(- alpha);

    //Forward eliminate
    for (int i = 1; i < N; i++){
        //Normalize row i (i in u convention):
        b(i-1) /= d(i-1);
        u(i) /= d(i-1); //Note: row i in u = row i-1 in the matrix
        d(i-1) = 1.0;
        //Eliminate
        u(i+1) += u(i)*alpha;
        d(i) += b(i-1)*alpha;
    }
    //Normalize bottom row
    u(N) /= d(N-1);
    d(N-1) = 1.0;

    // Backward substitute
    for (int i = N; i > 1; i--){ // loop from i=N to i=2
        u(i-1) -= u(i)*b(i-2);
        // b(i-2) = 0.0 // This is never read, why bother...
      }
}


void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */
    rowvec u_temp_row;
    for (int t = 1; t < T; t++){
        u.row(t) = u.row(t-1); // .copy() ?
        u_temp_row = u.row(t);
        tridiag(alpha,u_temp_row,N); //Note: Passing a pointer to row t, which is modified in-place
        u.row(t) = u_temp_row;

      }
}

void crank_nicolson(double alpha, mat &u, int N, int T){
    /*
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    */
    rowvec u_temp_row;
    rowvec u_temp_row1;
    for (int t = 1; t < T; t++){
        u_temp_row = u.row(t);
        u_temp_row1 = u.row(t-1);
        forward_step(alpha/2,u_temp_row,u_temp_row1,N);
        tridiag(alpha/2,u_temp_row,N);
        u.row(t) = u_temp_row;

      }
}
void g(mat &u, int N){
    // Initial condition u(x,0) = g(x), x \in (0,1)


    double dx = 1.0/N;

    for (int i = 1; i < N; i++){
      u(0,i) = 0;
  }
}

void analytic(mat &u, int N, int T, vec x, double dt){
  double L = 1.0;
  for (int t = 1; t < T; t++ ){
    for (int i = 1; i < N; i++){
      u(t,i) = x(i)/L - (2.0/(pi))*sin(x(i)*pi/L)*exp(-pi*pi*t*dt/(L*L));
    }
  }
}
