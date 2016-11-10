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
  params["num_branch_samples"] = 1;
  params["num_threads"] = 1;
  params["seed"] = 20120;  

  // Create these just for initializing
  int num_patterns = 2;
  double pattern_count_array[num_patterns];
  for(int i = 0; i < num_patterns; i++) pattern_count_array[i] = 1;

  // Set molecular model of sequence evolution
  double Q[16];
  string mol_model_name = "tamuraNei";
  makeQmatrix(mol_model_name, params, Q);

  map<int, string> bases;
  bases[0] = "T";
  bases[1] = "C";
  bases[2] = "A";
  bases[3] = "G";

  // Make some sequences
  vector<string> sequences;
  for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){
	sequences.push_back(bases[i] + bases[j]);
	//sequences.push_back(bases[i]);
      }
  }
  /*
  for(int i = 0; i < 16; i++)
    cout << sequences.at(i) << ' ';  
  cout << endl;
  */
  // Make a substitution model
  substModel mod(4, Q, "TCAG", num_patterns, pattern_count_array, params["relax"]);

  //Set up parallel random number streams
  createStreams(params["num_threads"], params["seed"]);  
  
  tree t;
  vector<double> logliks;
  double branch_one = 0.05;
  double branch_two = 0.05;
  vector<node_index> meristems;

  // Set up a file to write likelihood estimates to
  ofstream L_out;
  L_out.open("L.out", std::ios_base::app);

  // Now try matrix of possible transition probabilities
  logliks.clear();
  double num_reps = 100;
  double num_estimates = 1000;
  double branch_time = 0.05;
  double total_lik;
  for(int i = 0; i < 16; i++){
      for(int j = 0; j < 16; j++){
	L_out << sequences.at(i) << ' ' << sequences.at(j) << ' ';
	for(int h = 0; h < num_estimates; h++){
	  total_lik = 0;	 
	  for(int k = 0; k < num_reps; k++){
	    t = tree(mod, 0);
	    meristems.push_back(t.insert_lineage(0));
	    t.attach_sequence(meristems.at(0), sequences, branch_one, i, params);
	    meristems.push_back(t.insert_lineage(0));
	    t.attach_sequence(meristems.at(1), sequences, branch_two, j, params);
	    total_lik += exp(t.loglik());
	  }
	  L_out << total_lik / num_reps << ' ';
	}
	L_out << endl;
	//logliks.push_back(log(total_lik / num_reps));
	//cout << sequences.at(i) << ' ' << sequences.at(j) << ' ' << log(total_lik / num_reps) << endl;
      }
  }

  L_out.close();    

  /*
  cout << "\t T \t C \t A \t G \n";
  for(int i = 0; i < 4; i++){
    cout << bases[i] << '\t';
    for(int j = 0; j < 4; j++){
      cout << setprecision(4) << logliks.at(i*4 + j) << '\t';
    }
   cout << endl;
  }
  

  // Again try matrix of possible transition probabilities,
  // This time with a branch instead of both lineages going to the
  // root
  logliks.clear();
  for(int i = 0; i < 4; i++){
    sequences.at(0) = bases[i];
      for(int j = 0; j < 4; j++){
	sequences.at(1) = bases[j];
	total_lik = 0;
	for(int k = 0; k < num_reps; k++){
	  t = tree(mod, 0, TRUE); 
	  t.insert_root(0);
	  t.add_meristem (0, branch_time);
	  t.attach_sequence(0, sequences, branch_one + branch_time, 
			    0, params);
	  t.attach_sequence(1, sequences, branch_two + branch_time, 
			    1, params);
	  total_lik += exp(t.loglik());
	}
	logliks.push_back(log(total_lik / num_reps));
      }
  }

  cout << "\t T \t C \t A \t G \n";
  for(int i = 0; i < 4; i++){
    cout << bases[i] << '\t';
    for(int j = 0; j < 4; j++){
      cout << setprecision(4) << logliks.at(i*4 + j) << '\t';
    }
    cout << endl;
  }

  // Write the Q matrix to a file
  ofstream Q_out;
  Q_out.open("Q.txt", std::ios_base::app);
  for(int i = 0; i < 4; i++){
    for(int j = 0; j< 4; j++){
      Q_out << Q[i*4 + j] << '\t';
    }
    Q_out << endl;
  }

  // Write the likelihood matrix to a file
  ofstream L_out;
  L_out.open("L.txt", std::ios_base::app);
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      L_out << setprecision(4) << logliks.at(i*4 + j) << '\t';
    }
    L_out << endl;
  }
  L_out << endl;
  L_out.close();  
  */
  //Clean up
  destroyStreams();  
}
