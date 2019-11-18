#include "Functions.h"
#include "omp.h"


// Declare the function
void Probability(double Energy, vec &Energies, vec &counter);

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues, int &Nconfigs, bool randomconfig, vec &Energies, vec &counter)
{

  // Initialize the total number of accepted configurations
  Nconfigs = 0;
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment, randomconfig);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17);
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo cycles
  //#pragma omp parallel for
  for (int cycles = 1; cycles <= MCcycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int deltaE =  2*SpinMatrix(ix,iy)*
	  (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	   SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	   SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	   SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));

     // Metropolis test
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {

    // flip one spin and accept new spin config
	  SpinMatrix(ix,iy) *= -1.0;

    // Updating number of accepted configurations
    Nconfigs++;

    // Update energy and magnetisation
	  MagneticMoment += (double) 2*SpinMatrix(ix,iy);
	  Energy += (double) deltaE;


      }
    }
	}


  // Probability counting

  for (int i = 0; i < 400; i++){
    Energies(i) = -800 + 4*i;
  }

  if (MCcycles >= 10000){ // Hard coded the equilibrium state
    Probability(Energy, Energies, counter);
  }


    // update expectation values  for local node
    ExpectationValues(0) += Energy;    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;
    ExpectationValues(3) += MagneticMoment*MagneticMoment;
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins



// function to initialise energy, spin matrix and magnetization for ordered/unordered spin
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment, bool randomconfig)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // setup spin matrix
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){

      if (randomconfig == false){
      // setup spin matrix and initial magnetization
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
      }
    else {
    // Random orientation
      if ( RandomNumberGenerator(gen) >= 0.5 ){
        SpinMatrix(x, y) = 1.0;
            }
            else{
      	SpinMatrix(x, y) = -1.0;
            }
          }
  }
}
  // setup initial magnetization
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){

      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }

  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
	(SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
	 SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise


void WriteResultstoFile(ofstream& ofile, int NSpins, int MCcycles, double temperature, vec ExpectationValues, int Nconfigs, bool randomconfig)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;


  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  //ofile << "| Temperature | Energy-Mean | Magnetization-Mean|    Cv    | Susceptibility |\n";




  ofile << "\n";
  ofile << setw(20) << setprecision(8) << MCcycles; // # Monte Carlo cycles (sweeps per lattice)
  ofile << setw(20) << setprecision(8) << E_ExpectationValues/NSpins/NSpins; // Mean energy
  ofile << setw(20) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins; // Mean magetization
  ofile << setw(20) << setprecision(8) << Nconfigs*norm/NSpins/NSpins; // # accepted configurations
  ofile << setw(20) << setprecision(8) << Evariance/temperature/temperature; // Specific heat Cv
  ofile << setw(20) << setprecision(8) << Mvariance/temperature; // Susceptibility
  ofile << setw(20) << setprecision(8) << temperature;

} // end output function


void WriteResultsto4b(ofstream& ofile, int NSpins, int MCcycles, double temperature, vec ExpectationValues, int Nconfigs)
{ // divided by  number of MCcycles
  double norm = 1.0/((double) (MCcycles));
  
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;


  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues); // Energy Variance
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues); // Magnetization Variance

  double Estd = sqrt(Evariance*norm); // Standard deviation of the energy
  double Mstd = sqrt(Mvariance*norm); // Standard deviation of the magnetization

  ofile << setiosflags(ios::showpoint | ios::uppercase);

  ofile << "\n";
  ofile << setw(20) << setprecision(8) << MCcycles; // # Monte Carlo cycles (sweeps per lattice)
  ofile << setw(20) << setprecision(8) << E_ExpectationValues/NSpins/NSpins; // Mean energy
  ofile << setw(20) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins; // Mean magetization
  ofile << setw(20) << setprecision(8) << Evariance/temperature/temperature/NSpins/NSpins; // Specific heat Cv
  ofile << setw(15) << setprecision(8) << Mvariance/temperature/NSpins/NSpins; // Susceptibility
  ofile << setw(15) << setprecision(8) << Estd; // Standard deviation
  ofile << setw(15) << setprecision(8) << Mstd; // Standard deviation
} // end output function



void WriteConfigvsT(ofstream& ofile, int NSpins, int MCcycles, double temperature, int Nconfigs)
{
  double norm = 1.0/((double) (MCcycles));

  ofile << "\n";
  ofile << setw(20) << setprecision(8) << temperature;
  ofile << setw(20) << setprecision(8) << Nconfigs*norm/NSpins/NSpins; // # accepted configurations

} // end output function


void WriteResultstoFile2(ofstream& ofile, int NSpins, int MCcycles, double temperature, vec ExpectationValues, int Nconfigs)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;


  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  ofile << "\n";
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins; // Mean energy
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins; // Mean magetization
  ofile << setw(20) << setprecision(8) << Evariance/temperature/temperature; // Specific heat Cv
  ofile << setw(20) << setprecision(8) << Mvariance/temperature; // Susceptibility
} // end output function


// Function for counting the energy states => probability analysis for 4d

void Probability(double Energy, vec &Energies, vec &counter){
  double tol = 1E-10;
  for (int i = 0; i < 400; i++){
    if (fabs(Energy - Energies(i)) <= tol){
      counter(i) += 1;
    }
  }
}

// Writing a file for the probability output

void Writeprobabilities(ofstream& ofile, vec Energies, vec counter, int NSpins, int MCcycles, vec ExpectationValues)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;

  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues);

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  //ofile << "| Temperature | Energy-Mean | Magnetization-Mean|    Cv    | Susceptibility |\n";
  //ofile << "\n";
  for (int i = 0; i < 400; i++){
    ofile << setw(10) << setprecision(8) << Energies(i); // All the energies
    ofile << setw(15) << setprecision(8) << counter(i);
    ofile << "\n";
  }
  cout << "\n";
  cout << "Expectationvalue of the Energy = " << E_ExpectationValues << "\n"; // Mean Energy
  cout << "Variance of the Energy = " << Evariance << "\n"; // Variance
  cout << "Standard deviation of the Energy = " << sqrt(Evariance) << "\n"; // Standard Deviation

 // Probability distribution of the energy
} // end output function



// Write the results to the output file (expectation values)
void WriteT(ofstream& ofile, mat L, int NSpins, int MCcycles, vec T)
{

  // Divide by  number of Monte Carlo cycles
  double norm = 1.0/((double) (MCcycles));
  // Divide by number of spins
  double Norm = 1.0/(NSpins*NSpins);

  // Loop over temperature
  for (int i = 0; i < L.n_rows; i++){

    double E_ExpectationValues = L(i,0)*norm;
    double E2_ExpectationValues = L(i,1)*norm;
    double M_ExpectationValues = L(i,2)*norm;
    double M2_ExpectationValues = L(i,3)*norm;
    double Mabs_ExpectationValues = L(i,4)*norm;

    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
    double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;

    double Cv = Evariance/(T(i)*T(i));
    double Suscp = Mvariance/T(i);

    // Write to file
    ofile << setw(16) << setprecision(8) << T(i);
    ofile << setw(16) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
    ofile << setw(16) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins;
    ofile << setw(16) << setprecision(8) << Cv;
    ofile << setw(16) << setprecision(8) << Suscp;
    ofile << "\n";
  }

}
