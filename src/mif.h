// -*- mode: C++; -*-
#ifndef _MIF_H_
#define _MIF_H_

#include "lmatrix.h"
#include "substmodel.h"
#include "node.h"
#include "tree.h"
#include "usermodel.h"
#include "perturb.h"
#include <vector>
#include <list>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include "type_defs.h"


//typedef std::map<std::string, std::double>::iterator map_iterator;

class Mif {

private:
	
  double *condLoglik;
  vector<Usermodel> particles;

  class resampler {
    int n;
    double maxLoglik;
    double du;
    vector<int> tab;

  public:
    resampler (int np, double *w) : 
      n(np), tab(np) {
      
      //Scale log likelihoods by maximum log likelihood, transform to likelihoods
      maxLoglik = w[0];
      for (int j = 1; j < n; j++) 
        maxLoglik = (w[j] > maxLoglik) ? w[j] : maxLoglik;
      for (int j = 0; j < n; j++) 
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

    void resample (vector<Usermodel> &x, vector< map<string,double> > & y) {
      for (int i = 0, j = 0; i < n; i++) {
       if (tab[i] == 0) {
        while (tab[j] <= 1) j++;
         x[i] = x[j];
	 y[i] = y[j];
         tab[j]--;
        }
      }
    }
  };

public:
	
  //constructor
  Mif(substModel model, 
      string molModelName,
      int np, 
      int iterations,           
      int nNested, 
      int nThreads,
      map<string, double> & params,     // parameters
      map<string, double> & rw_sd,      // random walk standard deviations for each parameter
      map<string, double> & ivp_rw_sd,  // random walk standard deviations for IVPs
      vector<double> & times,           // times of observations 
      vector<string> & seqs,
      vector<double> & patternCount,
      string resDir,
      bool save_trees = FALSE,
      bool save_weights = FALSE)
    {
      //Extract the names of parameters to estimate
      vector<string> rw_params;
      rw_params.reserve(rw_sd.size() + ivp_rw_sd.size());
      for(map<string,double>::iterator it = rw_sd.begin(); it != rw_sd.end(); ++it) {
	rw_params.push_back(it->first);
      }
      for(map<string,double>::iterator it = ivp_rw_sd.begin(); it != ivp_rw_sd.end(); ++it) {
	rw_params.push_back(it->first);
      }

      //Weights and conditional log likelihood arrays
      double *weights = new double[np];
      //double *nestedWeights = new double[nThreads*nNested];
      int nt = times.size();
      condLoglik = new double[nt];
      
      //Initialize vector of parameters
      //NOTE: parameters are not perturbed in this step
      //AAK: could rewrite using matrices instead of maps
      vector< map<string, double> > param_vec; // parameters for each particle
      param_vec.reserve(np);
      for(int j = 0; j < np; j++)
	param_vec.push_back(params);

      // Reserve space for base number of particles
      // AAK: not all elements in seqs are sequences
      int numGnodes = 2*seqs.size() - 1;
      particles.reserve(np);

      // Build a perturb object for perturbing parameters
      Perturb perturber(rw_sd, ivp_rw_sd);

      // Perturb initial value parameters
      for (int j = 0; j < np; j++)
	perturber.perturb_initial_value_parameters(param_vec.at(j), 
						   param_vec.at(j)["cooling_factor"]);
      // Initialize particles
      for (int j = 0; j < np; j++)
	particles.push_back(Usermodel(model, param_vec.at(j), numGnodes));

      // Set up files to write results as the algorithm proceeds
      string resFileName = resDir + "res.txt";
      const char * resFile = resFileName.c_str();
      init_resfile(resFile, param_vec.at(0));
      string weightsFileName = resDir + "weights.txt";
      const char * weightsFile = weightsFileName.c_str();
      if(save_weights) initWeightsFile(weightsFile);

      // For timing
      time_t start, end;
      double iterTime;

      // Put the pattern count of each site pattern into an array
      // for lower level C codes to use
      // CHECK -- could rewrite to access patternCount.data (pointer)
      int numPatterns = patternCount.size();
      double patternCountArray[numPatterns];
      for(int i = 0; i < numPatterns; i++) patternCountArray[i] = patternCount.at(i);

      // RUN ITERATED FILTERING

      for (int k = 0; k < iterations; k++) { // Begin MIF loop
	
	// For timing
	start = time(0);

	// Initialize particles
	// AAK: is the order of perturbation ok for initial value parameters?
	if(k > 0) {
	  // Perturb initial value parameters
	  for (int j = 0; j < np; j++)
	    perturber.perturb_initial_value_parameters(param_vec.at(j), 
						       pow(param_vec.at(j)["cooling_factor"], k));

	  // Initialize particles
	  for (int j = 0; j < np; j++) {
	    particles.at(j) = Usermodel(model, param_vec.at(j), numGnodes);
	  }
	}

	for (int i = 0; i < nt; i++) { // Begin data loop

#pragma omp parallel for num_threads(nThreads) private(numPatterns)
	  for(int j = 0; j < np; j++) { // Begin particle loop

	    // Perturb parameters to be estimated		
	    perturber.perturb_parameters(param_vec.at(j), 
					 pow(param_vec.at(j)["cooling_factor"], k));

	    // Build a separate molecular model for each particle	    
	    double Q[16];
	    makeQmatrix(molModelName,param_vec.at(j),Q);
	    substModel newModel(4, Q, "TCAG", numPatterns, 
				patternCountArray, params["relax"]);

	    // Replace the old molecular model with the new
	    particles.at(j).set_mol_model(newModel);
	    
	    //Perform a nested resample of each particle
	    if (i == 0) { 
	      weights[j] = nestedResample(particles.at(j),
					  nNested,
					  param_vec.at(j), 
					  i,
					  param_vec.at(j)["start_time"],
					  times[i],
					  seqs,
					  newModel);	
	    } else {
	      weights[j] = nestedResample(particles.at(j),
					  nNested,
					  param_vec.at(j),
					  i,
					  times[i-1],
					  times[i],
					  seqs,
					  newModel);
	    }

	  } // End particle loop
	  
	  // Save weights to file on last iteration of MIF
	  if(save_weights && k == (iterations - 1)) 
	    updateWeightsFile(weightsFile, k, i, np, weights);

	  // Save genetic and infection trees to file if asked to
	  // Do so only for the last iteration of MIF
	  if(save_trees && (k == (iterations - 1)) && (i == (nt - 1))){
	    string temp = resDir + "gtrees.txt";
	    const char * gtree_file_name = temp.c_str();
	    for(int i = 0; i < np; i++) 
	      particles.at(i).save_gtree(gtree_file_name);
	    temp = resDir + "itrees.txt";
	    const char * itree_file_name = temp.c_str();
	    for(int i = 0; i < np; i++) 
	      particles.at(i).save_transmission_tree(itree_file_name);
	  }

	  //Systematic resampling
	  resampler R(np, weights);
	  R.resample(particles, param_vec);
	  condLoglik[i] = R.loglik();
	} // end data loop
	
	//For timing
	end = time(0);
	iterTime = end - start;

	//Write latest estimates to results file
	update_resfile(resFile, k, np, nt, rw_params,
		       param_vec, condLoglik, iterTime);

      } // end MIF loop
      
      //clean up
      delete[] weights;
    }
      
 
  //Write header of results file
  void init_resfile(const char * resFile, map<string,double> & params) {
    //open file
    ofstream out;
    out.open(resFile, std::ios_base::app);
    //Write header
    out << "iter ";
    map<string, double>::iterator iterator;
    for(iterator = params.begin(); iterator != params.end(); iterator++) {
      out << iterator->first << ' ';
    }
    out << "loglik runtime" << endl;
  }

  //Write results from latest MIF iteration  
  void update_resfile(const char * res_file, int iter, int np, int nt,
		      vector<string> rw_params,
		      vector< map<string,double> > & param_vec,
		      double * condloglik, double runtime){
    //open file
    ofstream out;
    out.open(res_file, std::ios_base::app);
    //Write results from this iteration
    out << iter << ' ';
    //Store parameter estimates from this MIF iteration
    double par_sum, par_avg;
    map<string, double>::iterator iterator;
    for(iterator = param_vec.at(0).begin(); 
	iterator != param_vec.at(0).end(); iterator++) {
      //Extract parameter name
      string par_name = iterator->first;
      if(std::find(rw_params.begin(), rw_params.end(), par_name) != rw_params.end()){
	//If the parameter is not fixed, average across the particles
	par_sum = 0;
	for (int j = 0; j < np; j++) { 
	  par_sum += param_vec.at(j)[par_name];
	}
	par_avg = par_sum / ( (double) np);
	out << par_avg << ' ';
      } else {
	//Otherwise, simply write the value of the fixed parameter
	out << param_vec.at(0)[par_name] << ' ';
      }
    }
    //Store log likelihood from this MIF iteration
    double loglik = 0;
    for(int i = 0; i < nt; i++) 
      loglik += condloglik[i];
    out << loglik << ' ' << runtime << endl;
  }

  //Set up file to store weights
  void initWeightsFile(const char * fileName){
    //open file
    ofstream out;
    out.open(fileName, std::ios_base::app);
    //Write header
    out << "mifIter dataPoint particle weight" << endl;
  }

  //Function to append this 
  void updateWeightsFile(const char * fileName, int mifIter, int dataPoint, int np, double * weights) {
    //open file
    ofstream out;
    out.open(fileName, std::ios_base::app);
    //Write weights
    for(int i = 0; i < np; i++) out << mifIter << ' ' << dataPoint << ' ' << i << ' ' << weights[i] << endl;
  }
      
  //Nested proposal and resample function
  double nestedResample(Usermodel & particle, const int & nExpand,
			map<string,double> & params, const int & sampleNum, 
			const double & timeStart, const double & timeEnd, vector<string> & seqs,
			substModel model) {

    //Allocate space for weights
    double nWeights[nExpand];

    // Reserve space for nested particles
    // AAK: reserving too many nodes and loop not needed if no sequence
    int numGnodes = 2*seqs.size() - 1;
    vector<Usermodel> nestedParticles;
    nestedParticles.reserve(nExpand);
    // Simulate to next sample
    particle.rprocess(params, timeStart, timeEnd);
    for (int j = 0; j < nExpand; j++)
      nestedParticles.push_back(particle);

    for(int i = 0; i < nExpand; i++){
      //Compute conditional log likelihood
      nWeights[i] = nestedParticles.at(i).dmeasure(params, seqs, sampleNum, timeEnd);
    }
    
    //Scale log likelihoods by maximum log likelihood, transform to likelihoods
    double maxLoglik = R_NegInf;
    for (int j = 0; j < nExpand; j++) {
      maxLoglik = (nWeights[j] > maxLoglik) ? nWeights[j] : maxLoglik;
    } 
    for (int j = 0; j < nExpand; j++) {
      nWeights[j] = exp(nWeights[j]-maxLoglik);
    }
    
     // AAK: check this
    if(!(R_FINITE(nWeights[nExpand -1]))) { // This is a check for NaN
       return R_NegInf; 
    }
    
    //Weighted sample of one particle; below selects the index of the sampled particle
    for (int j = 1; j < nExpand; j++) 
      nWeights[j] += nWeights[j-1];  //running sum
    double u = runif(0, nWeights[nExpand-1]);
    int sampledIndex = 0;
    while(u > nWeights[sampledIndex]) sampledIndex++;
    
    //Write sampled particle to particle vector
    particle = nestedParticles.at(sampledIndex);
    
    //Return particle weight
    double avgLik = nWeights[nExpand-1] / double(nExpand); //average likelihood   
    if (!(R_FINITE(log(avgLik) + maxLoglik))) {
       return R_NegInf;
    }
 
    return log(avgLik) + maxLoglik; //rescale to recover correct log liklihood
  }
  
  //destructor
  ~Mif (void) {
    if (condLoglik != NULL) 
      delete[] condLoglik;
  }

  //Extract conditional likelihood from ith sampling event
  double get_conditionalLoglik (int i) {
    return condLoglik[i];
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
};

#endif
