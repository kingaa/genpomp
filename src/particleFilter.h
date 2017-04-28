
// -*- mode: C++; -*-
#ifndef _PARTICLEFILTER_H_
#define _PARTICLEFILTER_H_

#include <algorithm>
#include "io.h"
#include "lmatrix.h"
#include "substmodel.h"
#include "node.h"
#include "tree.h"
#include "usermodel.h"
#include <vector>
#include <list>
#include <omp.h>
#include <fstream>
#include <iostream>
#include "type_defs.h"
#include "gsl_rng.h"
#include "gsl_randist.h"
#include "baseFilter.h"
//#include "mkl.h"

using namespace std;

class Particlefilter: public Basefilter {

private:
	
  double *condloglik;
  double *cond_hazard; // Components of the diagnosis likelihood
  double *cond_prob_no_diagnosis; // Components of the diagnosis likelihood
  vector<double> ess;
  vector<Usermodel> particles;
  vector< vector<int> > ancestors;
  vector< vector<double> > w;
  int nfail;
  	
public:
	
  // Constructor
  Particlefilter(gsl_rng * rngptr, 
		 substModel model, map<string,double> & params, 
		 vector<double> & times, vector<string> & seqs, 
		 string resdir, bool save_internals) {

    //mkl_disable_fast_mm();

    // Extract filtering parameters
    int np = params["np"];
    int num_nested = params["num_nested"];
    int num_threads = params["num_threads"];
    int nt = times.size();

    // Arrays for weights and conditional likelihoods
    double *weights = new double[np];
    double *hazards = new double[np];
    double *diagnosis_probs = new double[np];
      
    condloglik = new double[nt];
    cond_hazard = new double[nt]; // Components of the diagnosis likelihood
    cond_prob_no_diagnosis = new double[nt]; // Components of the diagnosis likelihood
    
    // Counter of number of filtering failures
    nfail = 0;

    /*
    //Reserve space for effective sample size and ancestor indices
    ess.reserve(nt);
    ancestors.reserve(nt);
    w.reserve(nt);
    */

    // Reserve space for base number of particles
    particles.reserve(np);
    for (int j = 0; j < np; j++)
      particles.push_back(Usermodel(model, params)); 

    // Write header if saving counts
    string states_file_name = resdir + "counts.txt";
    const char * states_file = states_file_name.c_str();
    if(save_internals) particles.at(0).write_states_header(states_file);


    // Run the particle filter
    cout << "Filtering with " << np << " particles over " << nt << " data points" << endl;    
    for (int i = 0; i < nt; i++) { 


      // If requested, write states to file
      if(save_internals){
	for(int k = 0; k < np; k++) particles.at(k).write_states(i, k, states_file);
      }

      std::vector<double> weights_vec;

      //mkl_cbwr_set(MKL_CBWR_AVX2);

      // BEGIN PARALLEL LOOP
#pragma omp parallel for num_threads(num_threads) schedule(static) private(weights_vec)
      for(int j = 0; j < np; j++) {

	if (i == 0) { 
	  weights_vec = nested_resample(&rngptr[j + 1], 
					particles.at(j), 
					num_nested, 
					params, 
					i, 
					params["start_time"], 
					times[i], 
					seqs, 
					model);	
	} else {
	  weights_vec = nested_resample(&rngptr[j + 1], 
					particles.at(j), 
					num_nested, 
					params, 
					i, 
					times[i - 1], 
					times[i], 
					seqs, 
					model);
	}
	hazards[j] = weights_vec.at(0); // these are on the log scale
	diagnosis_probs[j] = weights_vec.at(1); // these are on the log scale
	weights[j] = weights_vec.at(2); // cll
	// Ensure no memory leaks from mkl
	//mkl_free_buffers();
      }           
      // END PARALLEL LOOP      

      // Count failures
      for(int q = 0; q < np; q++){
	if(!R_FINITE(weights[q])) nfail += 1;
      }

      //Systematic resampling
      resampler R(&rngptr[0], np, weights);
      R.resample(particles);
      condloglik[i] = R.loglik();

      // Save averages of components of the diagnosis likelihood
      cond_hazard[i] = log_mean_exp(hazards, np); //NOTE: hazards is now modified
      cond_prob_no_diagnosis[i] = log_mean_exp(diagnosis_probs, np); //NOTE: diagnosis_probs is now modified
    }

    delete[] weights;
    delete[] hazards;
    delete[] diagnosis_probs;

    // Save conditional log likelihoods
    string res_file_name = resdir + "logliks.txt";
    const char * res_file = res_file_name.c_str();
    ofstream filestream;
    filestream.open(res_file);
    filestream << "hazard p_no_diagnosis condloglik" << endl;
    for(int j = 0; j < nt; j++) {
      filestream << cond_hazard[j] << ' ' << cond_prob_no_diagnosis[j] << ' ' << condloglik[j] << endl;
      }
    filestream.close();

    // Report number of filtering failures
    cout << nfail << " filtering failures" << endl;
  } 


  // Function to compute log mean exp on an array of numbers
  // The array is assumed to be on the log scale
  // NOTE: this function modifies the input
  double log_mean_exp(double *w, int n){
    // Compute the maximum of the array
    double max = w[0];
    for (int j = 1; j < n; j++) 
      max = (w[j] > max) ? w[j] : max;

    // Shift each element by the max and exponentiate
    for (int j = 0; j < n; j++) 
      w[j] = exp(w[j] - max);	

    // Compute the sum
    for (int j = 1; j < n; j++) 
      w[j] += w[j-1];

    // Average, take the log, then shift back by the max
    return log(w[n-1] / double(n)) + max;
  }

  //destructor
  ~Particlefilter (void) {
    if (condloglik != NULL) 
      delete[] condloglik;
    if (cond_hazard != NULL) 
      delete[] cond_hazard;
    if (cond_prob_no_diagnosis != NULL)
      delete[] cond_prob_no_diagnosis;
  }

  //Extract contional likelihood from ith sampling event
  double get_conditionalLoglik (int i) {
    return condloglik[i];
  }

  //Extract effective sample size
  vector<double> get_ess() {
    return ess;
  }

  //Extract weights
  vector< vector<double> > get_weights() {
    return w;
  }
  
  //Extract ancestor of each resampling event
  vector< vector<int> > get_ancestors() {
    return ancestors;
  }

};

#endif
