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
#include "padded_map.h"

//#include "mkl.h"

//extern int *thread_stream_map;
extern padded_map *thread_stream_map;

using namespace std;

class Particlefilter {

private:
	
  double *condloglik;
  double *cond_hazard; // Components of the diagnosis likelihood
  double *cond_prob_no_diagnosis; // Components of the diagnosis likelihood
  vector<double> ess;
  vector<Usermodel> particles;
  vector< vector<int> > ancestors;
  vector< vector<double> > w;
  int nfail;
  	
  class resampler {
    int n;
    double maxLoglik;
    double du;
    vector<int> tab;

  public:
    resampler (int np, double *w) : 
      n(np), tab(np) {
      
      // Copy the weights so they are not modified
      //double w[np];
      //for (int j = 1; j < np; j++) 
      //w[j] = weights[j];

      // Scale log likelihoods by maximum log likelihood, transform to likelihoods
      maxLoglik = w[0];
      for (int j = 1; j < np; j++) 
        maxLoglik = (w[j] > maxLoglik) ? w[j] : maxLoglik;

      // if(!R_FINITE(maxLoglik)) {
      // 	cout << "particle filter failure: all particles with weight -inf" << endl;
      // 	// EDIT: throw an exception here?
      // 	exit(1);
      // }

      for (int j = 0; j < np; j++) 
        w[j] = exp(w[j] - maxLoglik);	
      
      //Determine which indices are resampled
      for (int j = 1; j < n; j++) 
        w[j] += w[j-1];
      du = w[n-1] / double(n);
      double u = runif(-du,0);
    
      for (int i = 0, j = 0; i < n; i++) {
        u += du;
        while (u > w[j]) j++;
        tab[j]++;
      }
    }

    double loglik (void) {
      return log(du) + maxLoglik;
    }

    // Resample particles without tracking ancestors
    void resample (vector<Usermodel> &x) {
      for (int i = 0, j = 0; i < n; i++) {
	if (tab[i] == 0) {
	  while (tab[j] <= 1) j++;
	  x[i] = x[j];
	  tab[j]--;
        } 
      }
    }
    /*
    // Modifies ancestor vector
    void resample (vector<Usermodel> &x, vector< vector<int> > & ancestors) {
      vector<int> anc; // vector to hold indices of ancestors for each resampled particle
      anc.reserve(n);
      for (int i = 0, j = 0; i < n; i++) {
	if (tab[i] == 0) {
	  while (tab[j] <= 1) j++;
	  x[i] = x[j];
	  tab[j]--;
	  anc.push_back(j);
        } else {
	  anc.push_back(i);
	}
      }
      ancestors.push_back(anc);
    }
    */
  };

public:
	
  // Constructor
  Particlefilter(substModel model, map<string,double> & params, 
		 vector<double> & times, vector<string> & seqs, string resdir, bool save_internals) {

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
	
	// Tie the random number stream to the loop index
	thread_stream_map[omp_get_thread_num()].stream_index = j + 1;

	// std::vector<double> dummy;
	// dummy.push_back(runif(0,1));
	// dummy.push_back(runif(0,1));
	
	// if(runif(0,1) < 0.5){
	//   dummy.push_back(-1.0 /0.0);
	// } else {
	//   dummy.push_back(runif(0,1));
	// }
	// weights_vec = dummy;
	if (i == 0) { 
	  weights_vec = nested_resample(particles.at(j), num_nested, params, i, params["start_time"], times[i], seqs, model);	
	} else {
	  weights_vec = nested_resample(particles.at(j), num_nested, params, i, times[i-1], times[i], seqs, model);
	}
	hazards[j] = weights_vec.at(0); // these are on the log scale
	diagnosis_probs[j] = weights_vec.at(1); // these are on the log scale
	weights[j] = weights_vec.at(2); // cll
	// Ensure no memory leaks from mkl
	//mkl_free_buffers();
	
	// Reset rng stream to zeroth stream
	//thread_stream_map[omp_get_thread_num()] = 0;
      }           
      // END PARALLEL LOOP      

      for(int q = 0; q < num_threads; q++) thread_stream_map[q].stream_index = 0;

      // Count failures
      for(int q = 0; q < np; q++){
	if(!R_FINITE(weights[q])) nfail += 1;
      }
      
      // if(i == 6){
      // 	for(int q = 0; q < np; q++) cout << "weight " << q << " " <<  weights[q] << endl;
      // 	string tree_file_name = resdir + "trees.txt";
      // 	const char * tree_file = tree_file_name.c_str();
      // 	for(int q = 0; q < np; q++) particles.at(q).save_transmission_tree(tree_file);
      // 	string gtree_file_name = resdir + "gtrees.txt";
      // 	const char * gtree_file = gtree_file_name.c_str();
      // 	for(int q = 0; q < np; q++) particles.at(q).save_gtree(gtree_file);
      // }

      /*
      //Store the weights
      vector<double> tmp;
      tmp.reserve(np);
      for (int j = 0; j < np; j++) tmp.push_back(weights[j]);
      w.push_back(tmp);

      //Calculate and store the effective sample size for this obseravion
      double essDenom = 0;
      double denomJ;
      double sumWeights = 0;
      for (int j = 0; j < np; j++) if(isfinite(weights[j])) sumWeights += exp(weights[j]);
      for (int j = 0; j < np; j++) {
	denomJ = (exp(weights[j])/sumWeights)*(exp(weights[j])/sumWeights);
	if(isfinite(denomJ)) essDenom += denomJ;
      }
      ess.push_back(1/essDenom);
      */

      //Systematic resampling
      resampler R(np, weights);
      R.resample(particles);
      condloglik[i] = R.loglik();

      // if(R.loglik() > 1000){
      // 	stringstream ss;
      // 	ss << i;
      // 	string debug_file_name = resdir + "suspect_weights_data_point_" + ss.str() + ".txt";
      // 	const char * debug_file = debug_file_name.c_str();
      // 	ofstream filestream;
      // 	filestream.open(debug_file);
      // 	filestream << "index weight" << endl;
      // 	for(int q = 0; q < np; q++){
      // 	  filestream << q << ' ' <<  weights[q] << endl;
      // 	}
      // 	filestream.close();

      // 	string particle_file_name = resdir + "particle_" + ss.str() + ".txt";
      // 	const char * particle_file = particle_file_name.c_str();

      // 	for(int q = 0; q < np; q++){
      // 	  if( weights[q] > 1000){
      // 	    particles.at(q).print_particle(particle_file, params);
      // 	  }
      // 	}
      // }

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

  // Nested proposal and resample function
  std::vector<double> nested_resample(Usermodel & particle, const int & nExpand,
			map<string,double> & params, const int & sampleNum, 
			const double & timeStart, const double & timeEnd, vector<string> & seqs,
			substModel model) {

    // Simulate to next sample
    particle.rprocess(params,timeStart,timeEnd);
    
    // Extract diagnosis hazard and probability of no diagnosis 
      // (these will be the same for all nested particles)
    double log_hazard = particle.log_diagnosis_hazard(params);
    double log_prob_no_diagnosis = particle.log_prob_no_diagnosis(params);
    
    // Place above in results vector
    std::vector<double> results;
    results.push_back(log_hazard);
    results.push_back(log_prob_no_diagnosis);

    if(nExpand > 1){

      // Allocate space for weights
      double nWeights[nExpand];
      
      // Reserve space for nested particles
      //int numGnodes = 2*seqs.size() - 1;
      vector<Usermodel> nestedParticles;
      //nestedParticles.reserve(nExpand);
      for (int j = 0; j < nExpand; j++)
	nestedParticles.push_back(Usermodel(model, params));
            
      for(int i = 0; i < nExpand; i++){
	// Copy particle into buffer
	nestedParticles.at(i) = particle;
	// Compute conditional log likelihood
	nWeights[i] = nestedParticles.at(i).dmeasure(params,seqs,sampleNum,timeEnd);
      }
      
      // Find maximum of nested likelihoods
      double maxLoglik = nWeights[0];
      for (int j = 0; j < nExpand; j++) {
	maxLoglik = (nWeights[j] > maxLoglik) ? nWeights[j] : maxLoglik;
      } 
      
      //In the case of failures for all subsample particles, get out and return -INFINITY for loglik
      if(!(R_FINITE(maxLoglik))) { 
	particle = nestedParticles.at(0);
	results.push_back(NEG_INF);
	return results; 
      }

      //Scale log likelihoods by maximum log likelihood, transform to likelihoods
      for (int j = 0; j < nExpand; j++) {
	nWeights[j] = exp(nWeights[j]-maxLoglik);
      }
    
      //Weighted sample of one particle; below selects the index of the sampled particle
      for (int j = 1; j < nExpand; j++) 
	nWeights[j] += nWeights[j-1];  //running sum
      double u = runif(0,nWeights[nExpand-1]);
      int sampledIndex = 0;
      while(u > nWeights[sampledIndex]) sampledIndex++;
      
      //Write sampled particle to particle vector
      particle = nestedParticles.at(sampledIndex);

      //Return particle weight
      double avgLik = nWeights[nExpand-1] / double(nExpand); //average likelihood   
      
      if (!(R_FINITE(log(avgLik) + maxLoglik))) {
	results.push_back(NEG_INF);
	return results; 
      }
      
      results.push_back(log(avgLik) + maxLoglik); //rescale to recover correct log liklihood
      return results;
      
    } else {

      double weight = particle.dmeasure(params,seqs,sampleNum,timeEnd);
      results.push_back(weight); 
      return results;
      
    }
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

  /*

  //Nested proposal and resample function
  double nestedResample(Usermodel & particle, const int & nExpand,
			map<string,double> & params, const int & sampleNum, 
			const double & timeStart, const double & timeEnd, vector<string> & seqs,
			substModel model) {

    //Weights for nested particles
    double *nestedWeights = new double[nExpand];         
    //Reserve space for nested particles
    vector<Usermodel> nestedParticles;
    nestedParticles.reserve(nExpand);
    for (int j = 0; j < nExpand; j++)
      nestedParticles.push_back(Usermodel(model, params));

    //Simulate to next sample
    particle.rprocess(params,timeStart,timeEnd);
    for(int i = 0; i < nExpand; i++){
      //Copy particle into buffer
      nestedParticles.at(i) = particle;
      //Compute conditional log likelihood
      nestedWeights[i] = nestedParticles.at(i).dmeasure(params,seqs,sampleNum,timeEnd);
    }
    //Scale log likelihoods by maximum log likelihood, transform to likelihoods
    double maxLoglik = nestedWeights[0];
    for (int j = 1; j < nExpand; j++) 
      maxLoglik = (nestedWeights[j] > maxLoglik) ? nestedWeights[j] : maxLoglik;
    for (int j = 0; j < nExpand; j++) 
      nestedWeights[j] = exp(nestedWeights[j] - maxLoglik);	      
    
    //In the case of failures for all subsample particles, get out and return -INFINITY for loglik
    if(nestedWeights[nExpand -1] != nestedWeights[nExpand - 1]) { // This is a check for NaN
      delete[] nestedWeights;
      return -1.0/0.0; 
    }
    
    //Weighted sample of one particle; below selects the index of the sampled particle
    //RNGScope scope;
    for (int j = 1; j < nExpand; j++) 
      nestedWeights[j] += nestedWeights[j-1];  //running sum
    double u = runif(0,nestedWeights[nExpand-1]);
    int sampledIndex = 0;
    while(u > nestedWeights[sampledIndex]) sampledIndex++;
    
    //Write sampled particle to particle vector
    particle = nestedParticles.at(sampledIndex);
    
    //Return particle weight
    double avgLik = nestedWeights[nExpand-1] / double(nExpand); //average likelihood   
    delete[] nestedWeights;
    if( !isfinite(log(avgLik) + maxLoglik) ) return -1.0/0.0;
    return log(avgLik) + maxLoglik; //rescale to recover correct log liklihood
  }
  */
  
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
  /*
  //Extract tree from nth particle
  string get_tree (int n) {
    return particles.at(n).newick();
  }
  
  //Extract genetic tree from nth particle
  string get_gtree (int n) {
    return particles.at(n).gnewick();
  }
  */
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
