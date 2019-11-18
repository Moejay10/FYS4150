#ifndef FUNCTION_H
#define FUNCTION_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace  std;
using namespace arma;


// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&, bool);

// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &, int &, bool, vec &, vec &);

// The function to calculate the probability of the energy states.
void Probability(double, vec &, vec &);

// Prints to file the results of the calculations
// Task 4b
void WriteResultsto4b(ofstream&, int, int, double, vec, int);

// Task 4c
void WriteResultstoFile(ofstream&, int, int, double, vec, int, bool);

void WriteConfigvsT(ofstream&, int, int, double, int);

// Task 4d
void Writeprobabilities(ofstream&, vec, vec, int, int, vec);
// Task 4e
void WriteResultstoFile2(ofstream&, int, int, double, vec, int);

void WriteT(ofstream&, mat, int, int, vec);

#endif // ISING_H
