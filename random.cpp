#include <stdio.h>
#include <iostream>
#include <chrono>
#include "random.h"

//std::mt19937 my_rng {}; 
std::mt19937_64 my_rng {}; // Defines an 

std::uniform_real_distribution<double> my_unif_real_dist(0., 1.); //Define distribution
std::binomial_distribution<int> my_binomial_dist;

// Function to seed the random number generator from main file
// useful if you want the seed from a parameter file
// a negative value for seed gets you a random seed
// outputs the seed itself
int Seed(int seed)
{
  if (seed < 0) {
    long rseed=static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::cerr << "Randomizing random generator, seed is "<<rseed<<std::endl;
    my_rng.seed(rseed);
    return rseed;
  } else {
    std::cerr << "User-provided seed is "<<seed<<std::endl;
    my_rng.seed(seed);
    return seed;
  }
}
// This is the function to call if you want a random number in the interval [0,1)
double RANDOM(void)
{
  return my_unif_real_dist(my_rng);
}

int BINOMIAL(int N, double p)
{
  std::binomial_distribution<int> my_binomial_dist(N,p);
  return my_binomial_dist(my_rng);
}