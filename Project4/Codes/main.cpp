/*
   Program to solve the two-dimensional Ising model
   with zero external field and no parallelization
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
*/

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
// output file
ofstream ofile;



// Main program begins here

int main(int argc, char* argv[]) {

cout << "\n" << "Which Project Task do you want to run?: " << endl;
cout << "\n" << "Project Task A & B - Ising Model vs Analytical: " <<  "Write b " << endl;
cout << "\n" << "Project Task C - Equilibrium: " <<  "Write c " << endl;
cout << "\n" << "Project Task D - Probability Distribution: " <<  "Write d " << endl;
cout << "\n" << "Project Task E & F - Phase Transition & Critical Temperature: " <<  "Write e " << endl;


cout << "\n" << "Write here " << endl;
string Task;
cin >> Task;

if (Task == "b"){


// Analytical results for 2x2 lattice

double Z = 2*(exp(8) + exp(-8) + 6);
double E = 16*(exp(-8) - exp(8))/(Z);
double Mabs = 8*(exp(8) + 2)/(Z);
double Cv = 4*64*(4 + 6*exp(8) + 6*exp(-8))/(Z*Z);
double Xi = (32.0/(Z))*( (exp(8) + 1) - (2.0/(Z))* pow((exp(8) + 2),2) );

int NSpins = 2;

// Divide by number of spins
double norm = 1.0/NSpins/NSpins;

cout << "\n" << "Analytical values for all the expectation values divided by the number of spins \n";
cout << "Energy Mean per spin = " << E*norm << "\n";
cout << "Absolute Magnetization Mean per spin = " << Mabs*norm << "\n";
cout << "Specific heat per spin = " << Cv*norm << "\n";
cout << "Susceptibility per spin = " << Xi*norm << "\n";


cout << "\n" << "Project Task 4b: \n" << endl;

string file = "Lattice_2X2";
ofile.open(file);

ofile << "|   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility | Energy-STD | Magnetization-STD |\n";

int Nconfigs;
long int MCcycles;
double Temp = 1.0;

cout << "Read in the number of Monte Carlo cycles" << endl;
cin >> MCcycles;



vec Energies = zeros<mat>(400);
vec counter = zeros<mat>(400);


  // Start Monte Carlo sampling by looping over the selcted Temperatures
  //for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){

  for (int i = 10; i <= MCcycles; i *= 10){

  vec ExpectationValues = zeros<mat>(5);
  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(NSpins, i, Temp, ExpectationValues, Nconfigs, false, Energies, counter);

  WriteResultsto4b(ofile, NSpins, i, Temp, ExpectationValues, Nconfigs);
  }
 // }
  ofile.close();  // close output file
}


// Task 4c

if (Task == "c"){

  cout << "\n" << "Project Task 4c for ordered and unordered with spin L = 20: " << endl;

  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);

  string file = "Ordered";

  ofile.open(file);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "|   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean | # Accepted configurations |  Specific heat  | Susceptibility | Temperature |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N = 20;
  int Nconfigs;
  long int MC;
  double T;
  cout << "Read in the number of Monte Carlo cycles" << endl;
  cin >> MC;
  cout << "Read in the given value for the Temperature" << endl;
  cin >> T;


  //#pragma omp parallel for
  for (int i=1; i <= MC; i += 100){
    vec ExpectationValue = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, i, T, ExpectationValue, Nconfigs, false, Energies, counter);
    //
    WriteResultstoFile(ofile, N, i, T, ExpectationValue, Nconfigs, false);
  }
  ofile.close();  // close output file


  string file2 = "Unordered";

  ofile.open(file2);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| # Monte Carlo cycles  |  Energy-Mean   |  Magnetization-Mean  |  # Accepted configurations  |  Specific heat  |  Susceptibility   |  Temperature |\n";



  //#pragma omp parallel for
  for (int i=1; i <= MC; i += 100){
    vec ExpectationValue2 = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, i, T, ExpectationValue2, Nconfigs, true, Energies, counter);
    //
    WriteResultstoFile(ofile, N, i, T, ExpectationValue2, Nconfigs, true);
  }
  ofile.close();  // close output file

  // # Accepted configurations as a function of the Temperature

  string file3 = "Nconfig_vs_Temp";

  ofile.open(file3);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| Temperature | # Accepted configurations| \n";

  double InitialTemp = 1.0;
  double FinalTemp = 2.4;
  double TempStep = 0.2;

  for (double Temp = InitialTemp; Temp <= FinalTemp; Temp+=TempStep){
    vec ExpectationValue = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, 1000, Temp, ExpectationValue, Nconfigs, true, Energies, counter);

    WriteConfigvsT(ofile, N, 1000, Temp, Nconfigs);
   }
   ofile.close();  // close output file
}

// Task 4d

if (Task == "d"){

  cout << "\n" << "Project Task 4d for Probability, with spin L = 20 and T = 1.0: " << endl;

  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);

  string file = "Probability_1";

  ofile.open(file);

  ofile << "|  Energies | Energy counts |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N = 20;
  int Nconfigs; // # Accepted of configurations
  double T = 1.0; // Temperature

  long int MC;
  cout << "Read in the number of Monte Carlo cycles in times of 100" << endl;
  cin >> MC;

  int iterations;
  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false, Energies, counter);

    if (i == MC){
    Writeprobabilities(ofile, Energies, counter, N, iterations, ExpectationValue);
  }

  }
  ofile.close();  // close output file

  cout << "\n" << "Project Task 4d for Probability, with spin L = 20 and T = 2.4: " << endl;

  string file2 = "Probability_24";

  vec Energies2 = zeros<mat>(400);
  vec counter2 = zeros<mat>(400);

  T = 2.4; // Temperature
  ofile.open(file2);

  ofile << "|  Energies | Energy counts |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures

  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false, Energies2, counter2);

    if (i == MC){
    Writeprobabilities(ofile, Energies2, counter2, N, iterations, ExpectationValue);
  }

}
ofile.close();  // close output file

}

// Task 4c

if (Task == "e"){

cout << "Project Task 4e: " << endl;

vec Energies = zeros<mat>(400);
vec counter = zeros<mat>(400);

string file = "Temperature";

ofile.open(file);
ofile << "| Temperature | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";

// Start Monte Carlo sampling by looping over the selcted Temperatures
int N_start, N_step, N_final, n;
int Nconfigs;
long int MC;
double T_start, T_step, T_final, T;

// cout << "Read in the number of spins" << endl;
// cin >> N;
cout << "Read in the number of Monte Carlo cycles" << endl;
cin >> MC;


cout << "Read in the initial value for the Temperature" << endl;
cin >> T_start;
cout << "Read in the final value for the Temperature" << endl;
cin >> T_final;
cout << "Read in the number of integration points" << endl;
cin >> n;
T_step = (T_final - T_start)/((double) n);

cout << "T_step = " << T_step << endl;



// Declare a matrix which stores the expectation values for spins 40, 60, 80, 100
mat L_40 = zeros<mat>(n, 5);
mat L_60 = zeros<mat>(n, 5);
mat L_80 = zeros<mat>(n, 5);
mat L_100 = zeros<mat>(n, 5);

N_start = 40;
N_step = 20;
N_final = 100;

// Time the loop
double start = omp_get_wtime();
vec Tvalues = zeros<mat>(n);

for (int N = N_start; N <= N_final; N += N_step){



#pragma omp parallel for
// Start Monte Carlo sampling by looping over the selcted Temperatures
for (int i = 0; i < n; i++){
  vec ExpectationValue = zeros<mat>(5);


  T = T_start + T_step*i;
  Tvalues(i) = T;
  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(N, MC, T, ExpectationValue, Nconfigs, false, Energies, counter);
  //

  if (N == 40){
    L_40.row(i) = ExpectationValue.t();
  }

  else if (N == 60){
    L_60.row(i) = ExpectationValue.t();
  }

  else if (N == 80){
    L_80.row(i) = ExpectationValue.t();
  }

  else{
    L_100.row(i) = ExpectationValue.t();
  }

//#pragma omp ordered
  //WriteResultstoFile2(ofile, N, MC, T, ExpectationValue, Nconfigs);
}

}
double finish = omp_get_wtime();
double time_used = finish - start;
cout << "Time used [s]: " << time_used << endl;

WriteT(ofile, L_40, 40, MC, Tvalues);
//ofile <<"\n";
WriteT(ofile, L_60, 60, MC, Tvalues);
//ofile <<"\n";
WriteT(ofile, L_80, 80, MC, Tvalues);
//ofile <<"\n";
WriteT(ofile, L_100, 100, MC, Tvalues);
ofile.close();  // close output file




}

  return 0;

}
