// -*- mode: C++; -*-
#ifndef _BASE_FILTER_H_
#define _BASE_FILTER_H_

#include <algorithm>

class Basefilter {

private:

public:

  class resampler {
    int n;
    double maxLoglik;
    double du;
    vector<int> tab;

  public:

    resampler (gsl_rng * rngptr, int np, double *w) : 
      n(np), tab(np) {

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
      double u = gsl_runif(rngptr, -du, 0);
    
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
    // (for use with particle filter)
    void resample (vector<Usermodel> &x) {
      for (int i = 0, j = 0; i < n; i++) {
	if (tab[i] == 0) {
	  while (tab[j] <= 1) j++;
	  x[i] = x[j];
	  tab[j]--;
        } 
      }
    }

    // Overload resample member function
    // (for use with MIF)
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

  // Nested proposal and resample function
  std::vector<double> nested_resample(gsl_rng * rngptr, 
				      Usermodel & particle, 
				      const int & nExpand,
				      map<string,double> & params, 
				      const int & sampleNum, 
				      const double & timeStart, 
				      const double & timeEnd, vector<string> & seqs,
				      substModel model) {
  
    // Simulate to next sample
    particle.rprocess(rngptr, params, timeStart, timeEnd);
  
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
	nWeights[i] = nestedParticles.at(i).dmeasure(rngptr, params, seqs, sampleNum, timeEnd);
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
      double u = gsl_runif(rngptr, 0, nWeights[nExpand-1]);
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
    
      double weight = particle.dmeasure(rngptr, params, seqs, sampleNum, timeEnd);
      results.push_back(weight); 
      return results;
      
    }
  }

};

#endif




