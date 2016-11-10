// -*- mode: C++; -*-

/*

Constructs a simple tree, attaches sequences to that tree 
and computes the likelihood.

*/

#define MATHLIB_STANDALONE 
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <time.h>
#include <iomanip>
#include "lmatrix.h"
#include "userUnif.h"
#include "substmodel.h"
#include "tree.h"
#include "io.h"
#include "reduceSeqs.h"
#include "type_defs.h"

using namespace std;

int main ()
{

//80////////////////////////////////////////////////////////////////////////////////////

  // Constuct parameter map
  map<string, double> params;
  params["beta"] = 0.1;
  params["alphaR"] = 0.05;
  params["alphaY"] = 0.002;
  params["piT"] = 0.2;
  params["piC"] = 0.2;
  params["piA"] = 0.4;
  params["piG"] = 0.2;
  params["relax"] = 0;
  params["relax_branch"] = 0.1;
  params["stem_proportion"] = 0;
  params["fixed_stem"] = 0;
  params["num_branch_samples"] = 5;
  params["num_threads"] = 1;
  params["seed"] = 20999;  

  // Number of site patterns
  int numPatterns = 1;
  double patternCountArray[numPatterns];
  for(int i = 0; i < numPatterns; i++) patternCountArray[i] = 1;

  // Set molecular model of sequence evolution
  double Q[16];
  string molModelName = "tamuraNei";
  makeQmatrix(molModelName, params, Q);

  map<int, string> bases;
  bases[0] = "T";
  bases[1] = "C";
  bases[2] = "A";
  bases[3] = "G";
  // Print the Q matrix
  cout << "\n The Q Matrix: \n \n \t T \t C \t A \t G \n";
  for(int i = 0; i < 4; i++){
    cout << bases[i] << '\t';
    for(int j = 0; j< 4; j++){
      cout << Q[i*4 + j] << '\t';
    }
    cout << endl;
  }
  cout << endl;

  // Make a substitution model
  substModel mod(4, Q, "TCAG", numPatterns, patternCountArray, params["relax"]);

  // Make some sequences
  vector<string> sequences;
  sequences.push_back("T"); 
  sequences.push_back("C"); 
  sequences.push_back("G"); 

  //Set up parallel random number streams
  createStreams(params["num_threads"], params["seed"]);  
  
  tree t;
  vector<node_index> meristems;
  vector<double> logliks;
  vector<double> avg_logliks;
  int n = 1000;
  int k = 100;
  for(int j = 0; j < k; j++){  

    logliks.clear();
    // Make a tree with a single branch
    for(int i = 0; i < n; i++){
      t = tree(mod, 0);
      // Add three lineages
      meristems.push_back(t.insert_lineage(0));
      meristems.push_back(t.insert_lineage(0));
      // Branch a lineage
      meristems.push_back(t.add_meristem(meristems.at(0), 0.05));
      // Attach a sequence to each lineage
      t.attach_sequence(meristems.at(0), sequences, 0.15, 0, params);
      t.attach_sequence(meristems.at(1), sequences, 0.05, 1, params);
      t.attach_sequence(meristems.at(2), sequences, 0.15, 2, params);
      logliks.push_back(t.loglik());
    }
    /*
    cout << "Variability in one likelihood estimate when using " << 
      params["num_branch_samples"] << " gamma distributed branch length samples:" 
	 << endl << endl;
    for(int i = 0; i < 10; i++) cout << logliks.at(i) << endl;
    cout << endl;
    */

    double total = 0;
    for(int i = 0; i < logliks.size(); i++) total += exp(logliks.at(i));
    cout << "Estimate of log likelihood: " << log(total / logliks.size()) << 
      " (an average of " << n << " trees)" << endl;    
    avg_logliks.push_back(total / logliks.size());

  }

  // Print log likelihood
  double total = 0;
  for(int i = 0; i < avg_logliks.size(); i++) total += avg_logliks.at(i);
  cout << endl << "Average of estimates: " << log(total / avg_logliks.size()) <<
    endl;

  // Write likelihoods to file
  ofstream ll_out;
  ll_out.open("log_likelihoods.txt", std::ios_base::app);
  //ll_out  << log(total / avg_logliks.size()) << endl;
  for(int i = 0; i < avg_logliks.size(); i++) ll_out  << log(avg_logliks.at(i)) << endl;
  ll_out.close();
  
  //Clean up
  destroyStreams();  
}
