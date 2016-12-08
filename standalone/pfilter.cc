// -*- mode: C++; -*-

/*
Run a particle filter
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

#include "lmatrix.h"
#include "userUnif.h"
#include "substmodel.h"
#include "tree.h"
#include "abstractUsermodel.h"
#include "usermodel.h"
#include "particleFilter.h"
#include "io.h"
#include "reduceSeqs.h"
#include "type_defs.h"

using namespace std;

int main ( int argc, char *argv[] )
{
  if ( argc != 7) {
    cout << "usage: "<< argv[0] << " <paramsFile> <seqfile> <molModel> <resdir> <removeInvariantSites> <save_internals> \n";
    cout << "argc = " << argc << endl; 
  } else {

    // Set parameters
    map<string, double> params;
    readParams(params, argv[1]);

    // Read in sequences
    vector<double> times;
    vector<string> seqs;
    readSeqs(times, seqs, argv[2]);

    // Generate a reduced representation of the sequences if asked to
    bool removeInvariant;
    if(string(argv[5]) == "TRUE"){
      cout << endl <<  "Removing Invariant Sites" << endl << endl;
      removeInvariant = TRUE;
    } else if (string(argv[5]) == "FALSE"){
      cout << endl << "Not Removing Invariant Sites" << endl << endl;
      removeInvariant = FALSE;
    } else {
      throw(runtime_error("Please pass either \"TRUE\" or \"FALSE\" for <removeInvariantSites>"));
    }

    // Save internal workings of the filter if requested
    bool save_internals;
    if(string(argv[6]) == "TRUE"){
      cout << endl <<  "Saving internal workings of the particle filter" << endl << endl;
      save_internals = TRUE;
    } else if (string(argv[6]) == "FALSE"){
      cout << endl << "Running the particle filter without saving internal workings" << endl << endl;
      save_internals = FALSE;
    } else {
      throw(runtime_error("Please pass either \"TRUE\" or \"FALSE\" for <save_internals>"));
    }
    
    vector<double> patternCount = reduce_sequences(seqs, removeInvariant);
    int numPatterns = patternCount.size();

    double patternCountArray[numPatterns];
    for(int i = 0; i < numPatterns; i++) patternCountArray[i] = patternCount.at(i);

    // Set molecular model of sequence evolution
    double Q[16];
    string molModelName = argv[3];
    makeQmatrix(molModelName,params,Q);
    substModel mod(4, Q, "TCAG", numPatterns, patternCountArray, params["relax"]);

    // Set up parallel random number streams
    // One stream for the master and one for each iteration in the parallel loop
    //create_initial_rng_states(1 + params["np"], params["seed"]);  

    //Set up parallel random number generators
    gsl_rng * rngstreams = createStreams(1+ params["np"], params["seed"]);  

    //Run the particle filter
    Particlefilter pf(rngstreams, mod, params, times, seqs, argv[4], save_internals); 
    
    /*
    ofstream fileStream1;
    fileStream1.open(argv[4]);
    vector< vector<int> > ancs;
    ancs = pf.get_ancestors();
    for(int i = 0; i < ancs.at(0).size(); i++) { // np
      for(int j = 0; j < ancs.size() - 1; j++) { // nt - 1 
	fileStream1 << ancs.at(j).at(i) << ' ';
      }
      fileStream1 << endl;
    }
    fileStream1.close();

    //Save effective sample sizes
    ofstream fileStream2;
    fileStream2.open(argv[5]);
    vector<float_type> ess = pf.get_ess();
    for(int j = 0; j < ess.size(); j++) fileStream2 << ess.at(j) << endl;
    fileStream2.close();

    //Save Weights
    ofstream fileStream3;
    fileStream3.open(argv[6]);
    vector< vector<float_type> > w;
    w = pf.get_weights();
    for(int i = 0; i < w.at(0).size(); i++) { // np
      for(int j = 0; j < w.size(); j++) { // nt 
	fileStream3 << w.at(j).at(i) << ' ';
      }
      fileStream3 << endl;
    }
    fileStream3.close();

    //Save conditional log likelihoods
    ofstream fileStream4;
    fileStream4.open(argv[7]);
    for(int j = 0; j < seqs.size(); j++) 
      fileStream4 << pf.get_conditionalLoglik(j) << endl;
    fileStream4.close();
    */

    //Clean up
    destroyStreams();  
  }
}
