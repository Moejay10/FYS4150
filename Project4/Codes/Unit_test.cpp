#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "Functions.h"

TEST_CASE( "Checking the expectation values for a 2x2 lattice with T = 1.0, and MC-cycles = 10^6" ){

  // Define parameters for the test case
  int NSpins = 2;
  long int MCcycles = 1000000;
  int Nconfigs;
  double Temp = 1.0;


  // Divide by number of spins
  double norm = 1.0/NSpins/NSpins;

  // Analytical results for 2x2 lattice

  double Z = 2*(exp(8) + exp(-8) + 6);
  double E = 16*(exp(-8) - exp(8))/(Z);
  double Mabs = 8*(exp(8) + 2)/(Z);
  double Cv = 4*64*(4 + 6*exp(8) + 6*exp(-8))/(Z*Z);
  double Xi = (32.0/(Z))*( (exp(8) + 1) - (2.0/(Z))* pow((exp(8) + 2),2) );

  double E_Analytical = E*norm;
  double Mabs_Analytical = Mabs*norm;
  double Cv_Analytical = Cv*norm;
  double Xi_Analytical = Xi*norm;


  // Declare a matrix which stores the expectation values
  vec ExpectationValues = zeros<mat>(5);

  // Declare a vector which stores variables for equilibration and probability analysis
  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);

  // Run Monte Carlo with Metropolis sampling
  MetropolisSampling(NSpins, MCcycles, Temp, ExpectationValues, Nconfigs, false, Energies, counter);


  // divided by  number of MCcycles
  double Norm = 1.0/((double) (MCcycles));

  double E_ExpectationValues = ExpectationValues(0)*Norm;
  double E2_ExpectationValues = ExpectationValues(1)*Norm;
  double M_ExpectationValues = ExpectationValues(2)*Norm;
  double M2_ExpectationValues = ExpectationValues(3)*Norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*Norm;

  double E_Numerical = E_ExpectationValues*norm;
  double Mabs_Numerical = Mabs_ExpectationValues*norm;
  double Cv_Numerical = (E2_ExpectationValues - E_ExpectationValues*E_ExpectationValues)*norm/Temp/Temp;
  double Xi_Numerical = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)*norm/Temp;

  // Check the values
  REQUIRE( E_Numerical == Approx(E_Analytical).epsilon(0.001) );
  REQUIRE( Mabs_Numerical == Approx(Mabs_Analytical).epsilon(0.001) );
  REQUIRE( Cv_Numerical == Approx(Cv_Analytical).epsilon(0.1) );
  REQUIRE( Xi_Numerical == Approx(Xi_Analytical).epsilon(0.01) );
}
