#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <cstdlib>
#include "time.h"

using namespace std;
using namespace arma;
// Create an object for the output data file
ofstream ofile;
// Define functions for the source term and the closed-form solution
inline double f(double x) {return 100.0*exp(-10.0*x);}
inline double exact(double x) {return 1.0 - (1.0-exp(-10.0))*x - exp(-10.0*x);}


int main(int argc, char *argv[]){
  // Define a variable which stores the highest power of 10
  int exponent;
    // Read the highest power of 10 from the command line
    if( argc <= 1 ){
      cout << "Error: " << argv[0] << " reads max power of 10" << endl;
      exit(1);
    }
    else{
      exponent = atoi(argv[1]);
    }
    // Loop over the powers of 10
    for (int i = 1; i <= exponent; i++){
      // Give output data file the name data-i-.dat
      string fileout = "data";
      string argument = to_string(i);
      fileout.append(argument);
      fileout.append(".dat");

      // Define the number of grid points (excluding start and end points) and the step size
      int n = (int) pow(10.0,i);
      double h = 1.0/(n+1);
      double hh = h*h;

      // Define and set up an array for the x-values
      vec x(n+2);
      for (int i = 0; i <= n+1; i++) x(i) = i*h; 

      // GENERAL ALGORITHM ----------------------------------------------------------------------------------------------

      // Define arrays
      vec a_general(n); vec b_general(n); vec c_general(n); vec b_tilde_general(n); vec v_general(n+2);
      // Set up the endpoints of v
      v_general(0) = 0; v_general(n+1) = 0;
      // Setup of the decomposed matrix
      for (int i = 0; i <= n-1; i++){ 
	a_general(i) = -1.0;
	b_general(i) = 2.0;
	c_general(i) = -1.0;
	b_tilde_general(i) = hh*f(x(i+1));
      }

      // Run and time
      clock_t start, finish;
      start = clock();
      // Forward substitution
      for (int i = 1; i <= n-1; i++){
        double a_b = a_general(i-1)/b_general(i-1);
	b_general(i) -= a_b*c_general(i-1);
	b_tilde_general(i) -= a_b*b_tilde_general(i-1);
      }
      // Backward substitution
	v_general(n) = b_tilde_general(n-1)/b_general(n-1);
	for (int i = n-1; i >= 1; i--){
	  v_general(i) = (b_tilde_general(i-1) - c_general(i-1)*v_general(i+1))/b_general(i-1);
      }
      finish = clock();
      double time_general = (double) (finish - start)/(CLOCKS_PER_SEC);
      
      // SPECIAL ALGORITHM ---------------------------------------------------------------------------------------------

      // Define arrays
      vec b_special(n); vec b_tilde_special(n); vec v_special(n+2);
      // Set up the endpoints of v
      v_special(0) = 0; v_special(n+1) = 0;
      // Setup of the vectors b and b_tilde
      for (int i = 0; i <= n-1; i++){
	b_special(i) = (double) (i+2)/(i+1); 
	b_tilde_special(i) = hh*f(x(i+1));
      }

      // Run and time
      start = clock();
      // Forward substitution
      for (int i = 1; i <= n-1; i++) {
	b_tilde_special(i) += b_tilde_special(i-1)/b_special(i-1);
      }
      // Backward substitution
      v_special(n) = b_tilde_special(n-1)/b_special(n-1);
      for (int i = n-1; i >= 1; i--){
	v_special(i) = (b_tilde_special(i-1) + v_special(i+1))/b_special(i-1);
      }
      finish = clock();
      double time_special = (double) (finish - start)/(CLOCKS_PER_SEC);

      // LU-DECOMPOSITION -----------------------------------------------------------------------------------------------

      // Define matrix and arrays
      mat A = zeros<mat>(n,n);  vec b_tilde_lu(n); vec v_lu(n+2); 
      // Set up the endpoints of v
      v_lu(0) = 0; v_lu(n+1) = 0;
      // Setup of the matrix and b_tilde
      A(n-1,n-1) = 2.0; b_tilde_lu(n-1) = hh*f(x(n));
      for (int i = 0; i <= n-2; i++){
        A(i,i+1)  = -1.0;
        A(i,i)    = 2.0;
        A(i+1,i)  = -1.0;
	b_tilde_lu(i) = hh*f(x(i+1));
      }

      // Run and time
      start = clock();
      // Solve Av = b
      v_lu(span(1,n)) = solve(A,b_tilde_lu);
      finish = clock();
      double time_lu = (double) (finish - start)/(CLOCKS_PER_SEC);

      //------------------------------------------------------------------------------------------------------------------

      // Open output data file and write out results
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << "      x:             approx:        exact:         error:          "<< endl;
      for (int i = 0; i <= n+1 ;i++) {
	// Compute the relative error
	double epsilon = log10(fabs((v_special(i) - exact(x(i)))/exact(x(i))));
	ofile << setw(15) << setprecision(8) << x(i);
        ofile << setw(15) << setprecision(8) << v_special(i);
	ofile << setw(15) << setprecision(8) << exact(x(i));
	ofile << setw(15) << setprecision(8) << epsilon << endl;
      }
      
      // Write out log10 of the step size
      ofile << "log10(h):" << endl;
      ofile << setw(15) << setprecision(8) << log10(h) << endl;
      ofile.close();

      // Print the timing-results
      /*
      printf("n = 10^%d\n", i);
      printf("Time used by general algorithm: %.4f\n", time_general);
      printf("Time used by special algorithm: %.4f\n", time_special);
      printf("Time used by LU-decomposition algorithm: %.4f\n", time_lu);
      */
    }
    return 0;
}
