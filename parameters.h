#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace std;
class Parameter {
  public:
  //Parameter();
  //~Parameter(); //default destructor - unused

  long int random_seed =-1;
  int popsize = 50; //the total number of protocells
  int tau = 10; // the number of different ribozymes
  int numax = 8640; //the max size of a protocell, at this size it splits
  double epsilon=0.3;
  double constant1 = pow(tau,tau);
  double pdecay = 0.00000001; // probability that a ribozyme decays 
  string which_alpha_scheme="onehigher"; // onehigher,onelower,equal
                              // the different one is always
                              // the first

  int MaxTime=100000000; //default 100000000.;
  int time_save_data = popsize*10000; // how frequently we save data
  string output_filename = "data_pchrom.txt";
};

#endif

