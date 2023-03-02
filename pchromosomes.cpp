/*
 c++ version of the GARD model, trying to increase speed

TO DO:
- Try assign values to beta matrix from different distribution, maybe linear...
- Command line arguments
- read and write beta matrix - if beta matrix in command line arguments do not save to file again
- Ancestor's tree (instead of only the trace) should be implemented, 
because I see no coalescence within 1000s of generations

DONE:
- FUNCTION THAT PRINTS char *, ffs!!!  --- how about printf??? :P
- On line 376 there is a trick for optimisation that should make things 
A LOT FASTER !!! 
*/
#include <limits>
#include <iostream>
#include <fstream>      // std::ofstream, ifstream
#include <sstream>
#include <random>
#include <vector>
#include <math.h>       /* log,pow */
// #include <utility>      // std::pair, std::make_pair
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include "parameters.h"
#include "pchromosomes.h"
#include "random.h"

using namespace std;
extern Parameter par;
// CREATES POPULATION OF PROTOCELLS
vector <PROTOCELL> poppr( par.popsize , PROTOCELL(par.tau,par.numax/2) ); //the population of protocells

void CellDivision(vector<PROTOCELL>& poppr, int source_pos, int target_pos);

int main(int argc, char* argv[])
{
  
  int gtime,Time =0;
  double sumR;
  int sumvol;
  
  //check if file already exists, if it does the program does not start
  // this is here to avoid overwriting result files
  if(ifstream(par.output_filename)){
    std::cerr <<endl<< "Output file: "<< par.output_filename<<" already exists" << endl;
    std::cerr << "Simulation will not start" << endl << endl;
    exit(1);
  }
  
  Seed(par.random_seed); //randomises the seed of the random number generator
  int popsize = par.popsize;
    
  std::cerr << "par.time_save_data = "<<par.time_save_data << '\n';
  
  cout<<"Everything initialised. Beginning simulation!"<<endl<<endl;
  
  std::cerr << "popsize = "<<par.popsize<< ", tau = "<< par.tau << '\n';
  std::cerr << "maxtime = "<<par.MaxTime << '\n';
  
  // MAIN LOOP
  for(Time=0;Time<=par.MaxTime;Time++){
    // ANALYSIS AND PRINTING
    if(Time%(par.time_save_data)==0){
      // PrintData(Time, poppr);
      std::cerr << "Time: "<< Time << '\n';
      SaveData(Time, poppr);
    }
    
    for(gtime=0;gtime<par.popsize;gtime++){
      // REACTIONS
      //re-sets sumR 
      sumR=0.;
      sumvol=0;
      for(int i=0;i<popsize;i++){
        sumR += poppr[i].R;
        sumvol +=poppr[i].vol;
      }
      if(sumR==0. || sumvol==0){
        cout<<"Global extinction at Time = "<<Time<<", the program will terminate." << endl;
        exit(0);
        break;
      }
      
      // Replicate one ribozyme
      double rn = ((double)popsize)*RANDOM(); //pick a random number between 0 and sum(maxR) = popsize
      // double NONE=10.;
      // double rn = (NONE+sumR)*RANDOM();
      if(rn > sumR) continue;
      double counterR=0.;
      
      for(int i=0;i<popsize;i++){
        counterR += poppr[i].R;
        if(counterR>rn){
          //then protocell i is the one that replicates
          // which ribozyme will be copied? Assuming they all have same rate, pick one randomly
          // std::cerr << "Replicate rybozime, vol before = " << poppr[i].vol << '\n';
          
          poppr[i].ReplicateRibozyme(); //sets also volume and R
          // std::cerr << "after = " << poppr[i].vol << '\n';
          
          sumvol++;
          //if this cell is big enough it splits
          if(poppr[i].vol >= par.numax){
            //pick a random cell, different from current one
            int j=popsize*RANDOM();
            while(j==i){j=popsize*RANDOM();}
            sumvol -= poppr[j].vol;
            // std::cerr << "Cell division, vol j before = " << poppr[j].vol << '\n';
            
            CellDivision(poppr,i,j); //cell i splits into pos j
            // std::cerr << "after = " << poppr[j].vol << '\n';
            
            sumvol += poppr[j].vol;
          }
          break; //if we found and duplicated the ribozyme -> we exit this loop
        }
      }
      
      //let ribozymes decay
      // how many ribozymes decay per generation? binomial distribution
      int nr_riboz_decay = BINOMIAL(sumvol, par.pdecay/(double)popsize);
      // std::cerr << "nr_riboz_decay = "<< nr_riboz_decay << '\n';
      int rn2,countervol;
      for(int k=0;k<nr_riboz_decay;k++){
        rn2 = sumvol*RANDOM(); //pick a random ribozyme
        //now we have to find it
        countervol=0;
        for(int i=0;i<popsize;i++){
          countervol+=poppr[i].vol;
          if(countervol>rn2){
            // then we found the protocell, but which ribozyme?
            rn2 += -countervol + poppr[i].vol;
            poppr[i].DecayRibozyme(rn2); //decays one ribozyme at position rn
                                        // also takes care of change in R and vol  
            break; //we killed a ribozyme -> go to the next
          }
        }
        sumvol --;
      }
    }
  }
  
  cout<<"End of program"<<endl;
  
}
// END OF MAIN FUNCTION
void PrintData(int Time, std::vector<PROTOCELL> &poppr){
  std::cout << '\n'<<'\n';
  std::cout << "Time = " <<Time << '\n';
  for(int i=0;i<par.popsize;i++){
    poppr[i].PrintVol();
    poppr[i].PrintChromosomes();
    poppr[i].PrintR();
    std::cout << '\n';
  }
}

void SaveData(int Time, std::vector<PROTOCELL> &poppr){
  ofstream myfile;
  myfile.open(par.output_filename, std::ios_base::app);
  for(int i=0;i<par.popsize;i++){
    myfile << Time << " "<< poppr[i].R << " "<< poppr[i].vol;
    for(auto chr: poppr[i].chromosomes){
      myfile << " "<< chr; 
    }
    myfile << endl;
  }
  myfile.close();
  
}

// This function does cell division, cell at pos i splits and random half of content goes to j
// also takes care of 
void CellDivision(vector<PROTOCELL>& poppr, int source_pos, int target_pos){
  poppr[target_pos].DeleteCellContent(); //this cell dies, i.e. its content is deleted
  for(int i =0; i< par.tau;i++){
    int howmany_riboz = BINOMIAL( poppr[source_pos].chromosomes[i], 0.5 ); //random sample of ribozymes
    poppr[source_pos].chromosomes[i] -= howmany_riboz;
    poppr[target_pos].chromosomes[i] = howmany_riboz;
  }
  poppr[source_pos].Volume();
  poppr[target_pos].Volume();
  poppr[source_pos].MetabolicRate();
  poppr[target_pos].MetabolicRate();
}

// This is a PROTOCELL CLASS METHOD
// it duplicates one ribozyme
// it increases cell volume by 1
// it sets metabolic rate R

void PROTOCELL::ReplicateRibozyme(void){
  double totrepl=0.;
  for(int i=0;i<tau;i++) totrepl+=alpha[i]*chromosomes[i];
  double rn = totrepl*RANDOM(); //pick a random ribozyme, scaled by alpha
  double counterChr=0;
  // we have to find which ribozyme type this is
  for(int i=0;i<par.tau;i++){
    counterChr += alpha[i]*chromosomes[i];
    if(counterChr>rn){
      //we got the right ribozyme type to replicate
      chromosomes[i] ++ ;
      break;
    }
  }
  IncreaseVolBy(1); //increases volume by 1
  MetabolicRate(); // re-sets the metabolic rate
   
}

void PROTOCELL::DecayRibozyme(int pos){
  int counterChr=0;
  for(int i=0;i<par.tau;i++){
    counterChr += chromosomes[i];
    if(counterChr>pos){
      //we got the right ribozyme type to kill
      // std::cerr << "before dec chr["<<i<<"] = "<<chromosomes[i] << '\n';
      chromosomes[i] -- ;
      // std::cerr << "after dec chr["<<i<<"] = "<<chromosomes[i] << '\n';
      break;
    }
  }
  IncreaseVolBy( -1 ); //increases volume by -1 -> i.e. decreases it
  MetabolicRate(); // re-sets the metabolic rate
  // int vol1 = vol;
  // std::cerr <<endl<< "Volume after ribo decay: "<< vol << '\n';
  // Volume();
  // int vol2 = vol;
  // std::cerr << "Volume after recalc: "<< vol << '\n';
  // 
  // MetabolicRate(); // re-sets the metabolic rate
  // if(vol1!=vol2){
  //   std::cerr << "vol1!=vol2" << '\n';
  //   PrintChromosomes();
  //   PrintR();
  //   PrintVol();
  //   std::cerr << "now let's re-do the calculations" << '\n';
  //   Volume();
  //   MetabolicRate();
  //   PrintR();
  //   PrintVol();
  //   exit(1);
  // }
}