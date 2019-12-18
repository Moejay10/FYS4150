#include "Functions.h"
#include "omp.h"


// Code for solving the 1+1 dimensional diffusion equation
// du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
// with with L = 1, u(x,0) = g(x), u(0,t) = 0, u(L,t) = 1

double pi = 4*atan(1); // This is the constant pi

void forward_step(double alpha, rowvec &u, rowvec &uPrev, int N, bool CN){
    /*
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    */
    if (CN == false){
      for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (1.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
  }
    if (CN == true){
      for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (2.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
  }
}

void forward_step_2dim(double alpha, mat &u, mat &uPrev, int N){
    /*
    2-dimensional version of the forward-euler step. Uniform mesh.
    */

    for (int i = 1; i < N+1; i++){
      for (int j = 1; j < N+1; j++){
        u(i,j) = uPrev(i,j) + alpha*(uPrev(i+1,j) + uPrev(i-1,j)
        +uPrev(i,j+1)  + uPrev(i,j-1)- 4*uPrev(i,j));
      }
    }
}
void forward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements the forward Euler sheme, results saved to
    matrix u
    */

    // Skip boundary elements
    rowvec u_temp_row;
    rowvec u_temp_row1;
    for (int t = 1; t < T; t++){

        u_temp_row = u.row(t);
        u_temp_row1 = u.row(t-1);

        forward_step(alpha,u_temp_row,u_temp_row1,N,false);
        u.row(t) = u_temp_row;

      }
}
void forward_euler_2dim(double alpha, cube &u, int N, int T){
    /*
    2-dimensional version of forward Euler scheme, uniform mesh, results saved to cube u
    */

    // Skip boundary elements
    mat u_temp_mat;
    mat u_temp_mat1;
    for (int t = 1; t < T; t++){

        u_temp_mat = u.slice(t);
        u_temp_mat1 = u.slice(t-1);
        //cout << u.n_slices <<u.n_rows<<u.n_cols<< endl;

        forward_step_2dim(alpha, u_temp_mat,u_temp_mat1,N);
        u.slice(t) = u_temp_mat;

      }
}





void tridiagSolver(rowvec &u, rowvec u_prev, double alpha, int N, bool CN) {
  /*
  * Thomas algorithm:
  * Solves matrix vector equation Au = b,
  * for A being a tridiagonal matrix with constant
  * elements diag on main diagonal and offdiag on the off diagonals.
  */
 double diag, offdiag;
 if (CN == false){
   diag = 1+2*alpha;
}

 if (CN == true){
  diag = 2+2*alpha;
}

 offdiag = -alpha;
 vec beta = zeros<vec>(N+1); beta(0) = diag;
 vec u_old = zeros<vec>(N+1); u_old(0) = u_prev(0);
 double btemp;


 // Forward substitution
 for(int i = 1; i < N+1; i++){
   btemp = offdiag/beta[i-1];
   beta(i) = diag - offdiag*btemp;
   u_old(i) = u_prev(i) - u_old(i-1)*btemp;
 }


 // Special case, boundary conditions
 u(0) = 0;
 u(N+1) = 1;

 // backward substitution
 for(int i = N; i > 0; i--){
   u(i) = (u_old(i) - offdiag*u(i+1))/beta(i);
  }

}





void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */
    rowvec u_temp;
    rowvec u_temp1;
    for (int t = 1; t < T; t++){
        u_temp = u.row(t);
        u_temp1 = u.row(t-1);
        tridiagSolver(u_temp, u_temp1, alpha, N, false); //Note: Passing a pointer to row t, which is modified in-place
        u_temp(0) = 0;
        u_temp(N+1) = 1;
        u.row(t) = u_temp;
      }
}

void crank_nicolson(double alpha, mat &u, int N, int T){
    /*
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    */
    rowvec u_temp;
    rowvec u_temp1;
    rowvec u_temp2;
    for (int t = 1; t < T; t++){
      u_temp = u.row(t-1);
      u_temp1 = u.row(t);
      forward_step(alpha, u_temp1, u_temp, N, true);

      u_temp1(0) = 0;
      u_temp1(N+1) = 1;

      u_temp2 = u.row(t);
      tridiagSolver(u_temp2, u_temp1, alpha, N, true);
      u_temp2(0) = 0;
      u_temp2(N+1) = 1;
      u.row(t) = u_temp2;
      }
}

void analytic(mat &u, int N, int T, int infty) {
  /*
   * Analytic solution for the one dimensional diffusion equation,
   * with boundary conditions:
   * u(L,t) = 1 and u(0,t) = 0, for all t.
   * u(x,0) = 0 for x < L.
   */
  double L = 1.0; // scale the rod such that x goes from 0 to L=1.
  double sum;
  vec x = linspace<vec>(0,1,N+2);
  vec t = linspace<vec>(0,1,T);

  // time loop
  for (int i = 0; i < T; i++){
    // position x loop
    for (int j = 0; j < N+2; j++){
      // calculate the transient solution.
      sum = 0;
      for (int n = 1; n < infty; n++){
        sum +=  ((2*pow(-1,n))/n/pi)*sin(n*pi*x(j)/L)*exp(-n*n*pi*pi*t(i)/L/L);
      }
      u(i,j) = x(j)/L + sum;
      //cout << "Sum in Analytical solution = " << sum << endl;
    }
  }
}


void analytic_2D(mat &u, int N, double t){
  /*
   * Analytic solution for the two dimensional diffusion equation at a
   * given time. Boundary conditions are all equal to zero. A "sine-paraboloid"
   * which is centered in the middle of the medium is the initial condition.
   * Boundary conditions are zeros.
   */
  double L = 1.0;
  vec x = linspace<vec>(0,1,N+2);
  vec y = linspace<vec>(0,1,N+2);


  for (int i = 0; i < N+2; i++){
    for (int j = 0; j < N+2; j++){
      u(i,j) = sin(pi*x(i))*sin(pi*y(j))*exp(-2*pi*pi*t);
    }
  }
}



// Function for setting up the iterative Jacobi solver
int JacobiSolver(mat &u, double dx, double dt, double abstol, int maxiter){
  ofstream ofile;

  string dxstring = to_string(dx);
  dxstring.erase ( dxstring.find_last_not_of('0') + 1, std::string::npos );
  string dtstring = to_string(dt);
  dtstring.erase ( dtstring.find_last_not_of('0') + 1, std::string::npos );
  string file = "2dim_Implicit:dx="+dxstring+"dt="+dtstring;
  ofile.open(file);
  ofile << u;

  mat u_prev = u;

// Constants, variables
  int N = int(1/dx);   // Number of integration points along x & y -axis (inner points only)
  int T = int(1/dt);   // Number of time steps till final time

  double alpha = dt/(dx*dx);
  double term = 1.0/(1.0 + 4*alpha);
  double term1 = alpha*term;
  double scale = (N+2)*(N+2);
  double Delta;
  double diff;
  int i,j;


// Time loop
// Parallelization using openmp
//#pragma omp parallel for
for (int t = 1; t < T; t++){
  int iterations = 0;
  diff=1;
  mat u_guess = ones<mat>(N+2,N+2);

  while (iterations < maxiter && diff > abstol){
    diff = 0;
    // Define parallel region
    // Loop over all inner elements, will converge towards solution
    for (j = 1; j < N+1; j++) {
      for (i = 1; i < N+1; i++) {
        // u_guess is the previous u, which also work for a random guess
        Delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
        u(i,j) = term1*Delta + term*u_prev(i,j);
        diff += fabs(u(i,j) - u_guess(i,j));
          }
        } // end of double for loop
    u_guess = u;
    diff /= scale;
    iterations++;
    }  // end iteration loop
  u_prev = u;
  ofile << u;
  } // end time loop
ofile.close();
}


void Lithosphere(int Case, double dx, double dt, double tol, int maxiter){
  /*
   1 Ga of the lithosphere in three zones:
   Zone 1: 0-20km,
   Zone 2: 20-40km,
   Zone 3: 40-120km.
   Case=1: No heat production
   Case=2: Heat production
   Case=3: Mantle and no decay
   Case=4: Mantle with decay
   */

   ofstream ofile;
   string file;


  double Qmantle = 0; // represents the constant Q term in the mantle
  if (Case == 1) {
    file = "No_Heat";
  } else if (Case == 2) {
    file = "Heat";
    Qmantle = 0.05;
  } else if (Case == 3) {
    file = "No_Decay";
    Qmantle = 0.55;
  } else if (Case == 4) {
    file = "Decay";
    Qmantle = 0.05;
  } else {
    cout << "Error, Case must be 1,2,3,4.\n";
    exit(1);
  }
  ofile.open(file);
  // initialization
  int Nx = int(1.25/dx); // 0-150 km wide
  int Ny = int(1/dx); // 0-120 km deep
  int T = int(1/dt); // number of time steps
  mat u = zeros<mat>(Nx+1,Ny+1);

// Applying the cases
  if (Case == 1){
    NoHeat(u,Nx,Ny);
  }
  else {
    Heat(u,Nx,Ny);
  }
  mat u_prev = u;

  // Constants and variables
  double Delta, diff, Q; // temporary constant that will be used
  double scale = (Nx+1)*(Ny+1); // No. of points in matrix
  double alpha = dt/(dx*dx);
  double rho = 3.51*1e3; // Density [kg/m^3]
  double cp = 1000; // Heat capacity [J/kg/deg Celsius]
  double k = 2.5; // Thermal conductivity [W/m/deg Celsius]
  double T_s = 1292; // Temperature scale - T(120) - T(0)
  double t_s = 3.1557e16; // Time scale
  double x_s = 120000; // Distance scale - xmax = 120 km
  double beta = t_s*k/(rho*cp*x_s*x_s);
  double Q_s = 1e-6*t_s*dt/(rho*cp*T_s);
  double term1 = beta*alpha;
  double term2 = 1/(1 + 4*term1);

  // Construction of Q_vec: stores constant heat production in zones:
  vec Q_vec = zeros<vec>(Ny+1);
  if (Case != 1){
    double Q1 = 1.4*Q_s;
    double Q2 = 0.35*Q_s;
    double Q3 = Qmantle*Q_s;
    for (int i = 0; i < 17; i++) Q_vec(i) = Q1; // 0-20 km depth
    for (int i = 17; i < 34; i++) Q_vec(i) = Q2; // 20-40 km depth
    for (int i = 34; i <= Ny; i++) Q_vec(i) = Q3; // 40-120 km depth
  }

  // EConstruction of Q_time: stores time dependet heat production (due to decay):
  vec Q_time = zeros<vec>(T+1);
  if (Case == 4){ // Q_time is simply 0 for cases 1-3
    Decay(Q_time, T, Q_s);
  }


  // Time loop
  for (int t = 1; t <= T; t++){
    int iter = 0;
    double diff=1;
    mat u_guess = u; // Random matrix, doens't matter what, will converge
    // Jacobi method loop
    while (iter < maxiter && diff > tol){
      diff = 0;
      for (int j = 1; j < Ny; j++){
          // Case separation
          Q = Q_vec(j); // constant Q for depth at index j
          if (j >= 34) Q += Q_time(t); // time dependent Q, only at depth j>=34
        for (int i = 1; i < Nx; i++){
          // Jacobi method algorithm, u should converge towards solution
          Delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
          u(i,j) = term2*(term1*Delta + Q + u_prev(i,j));
          diff += fabs(u(i,j) - u_guess(i,j));
        }
      }
      u_guess = u;
      diff /= scale;
      iter++;
    } // end Jacobi method loop
    u_prev = u;
    //cout << "timestep: " <<  t << " iterations: " << iter << endl;
    ofile << u;
  } // end time loop
  ofile.close();
}

void NoHeat(mat& u, int Nx, int Ny){
  // Creates boundary conditions for the case of Q = 0 everywhere
  for (int i = 0; i < Nx+1; i++){
    u(i,0) = 0;
    u(i,Ny) = 1;
  }
  // linear temperature from 0-1 (scaled)
  for (int i = 0; i < Nx+1; i += Nx){
    for (double j = 0; j < Ny+1; j++){
      u(i,j) = (double) j/Ny;
      u(i,j) = (double) j/Ny;
    }
  }
}

void Heat(mat& u, int Nx, int Ny){
  // Analytical solution to steady state Temp(depth) with heat production
  // Temperature is scaled [8-1300] -> [0-1292] -> [0,1]
  double z, func, limit;
  double a_1, a_2, a_3; // Coefficients for the functions
  double b_1, b_2, b_3; // Coefficients for the functions
  double c_1, c_2, c_3; // Coefficients for the functions

  a_1 = -0.28; a_2 = -0.07; a_3 = -0.01;
  b_1 = 23.66; b_2 = 15.26; b_3 = 10.46;
  c_1 = 0.0; c_2 = 86; c_3 = 180;

  limit = 17;

  for (double j = 0; j < limit; j++){
    z = j*1.2;
    func = (a_1*z*z + b_1*z + c_1)/1292.0;
    for (int i = 0; i <= Nx; i++){
    u(i,j) = func;
    }
  }
  for (double j = limit; j < 2*limit; j++){
    z = j*1.2;
    func = (a_2*z*z + b_2*z + c_2)/1292.0;
    for (int i = 0; i <= Nx; i++){
      u(i,j) = func;
    }
  }
  for (double j = 2*limit; j <= Ny; j++){
    z = j*1.2;
    func = (a_3*z*z + b_3*z + c_3)/1292.0;
    for (int i = 0; i <= Nx; i++){
      u(i,j) = func;
    }
  }
}

void Decay(vec& Q_time, int T, double Q_s){
  double Uranium, Thorium, Potassium, times;
  for (int t = 0; t <= T; t++){
    times = (double) t/T;
    Uranium = exp (-0.155*times);
    Thorium = exp (-0.0495*times);
    Potassium = exp (-0.555*times);
    Q_time(t) = Q_s*(0.2*Uranium + 0.2*Thorium + 0.1*Potassium);
  }
}
