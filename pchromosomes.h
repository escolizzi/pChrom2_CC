#ifndef PROTOCELL_H
#define PROTOCELL_H

#include <iostream>
#include <string>
#include <vector>

#include "parameters.h"
#include "random.h"
using namespace std;
extern Parameter par;

//class PROTOCELL, contains all the information about protocells
class PROTOCELL{
public:
  // variables and data
  int vol; //volume, i.e. total count of ribozymes
  double R; //metabolic rate of the cell
  int tau=par.tau;
  // choromosomes contains the ribozyme types, and their abundance
  vector<int> chromosomes;
  vector<double> alpha;

  // default creator, for an empty protocell
  // PROTOCELL()
  // {
  //   vol=0;
  //   R=0.;
  //   tau = par.tau;
  // }
  //creator if you want an equal number of ribozymes 
  // for each ribozyme class
  PROTOCELL(int tau, int init_vol)
  {
    tau = tau;
    for(int i=0;i<tau;i++){
      chromosomes.push_back(init_vol/tau);
      if(par.which_alpha_scheme=="onehigher"){
        if(i==0) alpha.push_back(1.0);
        else alpha.push_back(0.9);
      }else if(par.which_alpha_scheme=="onelower"){
        if(i==0) alpha.push_back(1.0);
        else alpha.push_back(0.9);
      }else if(par.which_alpha_scheme=="equal"){
        alpha.push_back(1.0);
      }else{
        cerr<<"Error: which_alpha_scheme = "<< par.which_alpha_scheme << " not recognised.";
        exit(1);
      }

    }
    Volume();
    MetabolicRate();
  }
  
  // destructor
  ~PROTOCELL(){}
    
  // calculates volume of a protocell
  void Volume(){
    vol=0;
    for(auto chr: chromosomes){
      vol+=chr;
    }
  }
  //calculates metabolic rate of protocell
  void MetabolicRate(){
    if(vol==0){
      R=0.;
    }else{
      R = 1.;
      for(auto chr: chromosomes){
        R = R* par.tau*chr/(double)vol;
      }
      R = (vol/(double)par.numax) * pow(R , par.epsilon);
    }
  }
  
  //picks a random ribozyme for replication
  //also increases volume and sets R
  void ReplicateRibozyme(void);
  void DecayRibozyme(int pos);
  
  // functions that act on protocells
  void IncreaseVolBy(int howmuch){
    vol += howmuch;
  }
  
  void DeleteCellContent(void){
    for(int i=0;i<tau;i++){
      chromosomes[i]=0;
    }
    vol=0;
    R=0.;
  }
  
  void PrintVol(void){
    cout<<"Vol = "<< vol<<endl;
  }
  
  void PrintChromosomes(void){
    for(int i=0;i<tau;i++){
      std::cout<<i<<","<<chromosomes[i]<<" ";
    }
    cout << endl;
  }
  void PrintR(void){
    cout<<"R = "<<R<<endl;
  }
};

void PrintData(int Time, std::vector<PROTOCELL> &poppr);
void SaveData(int Time, std::vector<PROTOCELL> &poppr);

#endif