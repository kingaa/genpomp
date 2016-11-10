// -*- mode: C++; -*-

/*
Constructs a simple tree, simulates sequences on that tree.
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

using namespace std;



int main ( int argc, char *argv[] )
{
  if ( argc != 9) {
    cout << "usage: "<< argv[0] << " <nsims> <piT> <piC> <piA> <piG> <alphaR> <alphaY> <beta> \n";
    cout << "argc = " << argc << endl; 
  } else {

//80////////////////////////////////////////////////////////////////////////////////////

    // Constuct parameter map
    map<string,double> params;
    params["beta"] = strtod(argv[8], NULL);;
    params["alphaR"] = strtod(argv[6], NULL);;
    params["alphaY"] = strtod(argv[7], NULL);;
    params["piT"] = strtod(argv[2], NULL);
    params["piC"] = strtod(argv[3], NULL);
    params["piA"] = strtod(argv[4], NULL);;
    params["piG"] = strtod(argv[5], NULL);;
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
      }
    }
    
    // Make a substitution model
    substModel mod(4, Q, "TCAG", num_patterns, pattern_count_array, params["relax"]);
    
    //Set up parallel random number streams
    createStreams(params["num_threads"], params["seed"]);  
    
    tree t;
    double branch_one = 0.05;
    double branch_two = 0.05;
    double branch_three = 0.05;
    vector<node_index> meristems;
    
    // Set up a file to write simulated sequences to 
    ofstream seqs_out;
    seqs_out.open("seqs.out", std::ios_base::app);
    
    // Simulate sequences
    double branch_length = 1;
    int nsims = atoi(argv[1]);
    for(int i = 0; i < nsims; i++){
      t = tree(mod, 0);
      // Add three lineages
      meristems.push_back(t.insert_lineage(0));
      meristems.push_back(t.insert_lineage(0));
      // Branch a lineage
      meristems.push_back(t.add_meristem(meristems.at(0), branch_length/2));
      // Simulate a sequence on each meristem
      seqs_out << t.simulate_sequence(meristems.at(0), 2, branch_length*1.5, params) << ' ';
      seqs_out << t.simulate_sequence(meristems.at(1), 2, branch_length/2, params) << ' ';
      seqs_out << t.simulate_sequence(meristems.at(2), 2, branch_length*1.5, params) << endl;
    }
    
    seqs_out.close();    
    destroyStreams();  
  }
}
