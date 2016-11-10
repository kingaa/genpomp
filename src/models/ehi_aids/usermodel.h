// -*- C++ -*-
//EHI = 'Early HIV Infection'
//This model has six compartments (3 stages of disease by 2 diagnosis states)
//This new formulation of the early HIV model does not assume constant incidence
//The likelihood function has two components: a genetic likelihood and a diagnosis likelihood

#ifndef _USERMODEL_H_
#define _USERMODEL_H_

#include <R.h>
#include <Rmath.h>
#include "abstractUsermodel.h"
#include "lmatrix.h"
#include "substmodel.h"
#include "node.h"
#include "tree.h"
#include <fstream>
#include <vector>
#include <list>
#include <stdexcept>
#include "type_defs.h"

using namespace std;
 
class Usermodel : public AbstractUsermodel
{

private:

  //vectors of meristem indices in each class of infected individuals
  vector<node_index> _undiagnosed0;	
  vector<node_index> _undiagnosed1;
  vector<node_index> _undiagnosed2;
  vector<node_index> _diagnosed0;
  vector<node_index> _diagnosed1;  
  vector<node_index> _diagnosed2;  
  
  //counts of individuals of each class
  vector<int> _I0;
  vector<int> _I1;
  vector<int> _I2;
  vector<int> _J0;
  vector<int> _J1;      
  vector<int> _J2;      

  //event times
  vector<double> _times;

  //sequences
  vector<string> _seqs;
  vector<double> _sequence_times;

  //accumulator variables to integrate counts of I0, I1, and I2 individuals
  double _I0_pop_integral;
  double _I1_pop_integral;
  double _I2_pop_integral;

  //Vector of states for plotting with ouch
  vector<string> _states;

public:

  //constructor
  Usermodel(substModel model, 
	    map<string,double> & params,
	    bool make_tree_table = FALSE,
	    int num_gnodes = 0) : AbstractUsermodel(model, 
						    params, 
						    make_tree_table,
						    num_gnodes)
  {
    for (int i = 0; i < params["numI0"]; i++) {
      _undiagnosed0.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("I0");
    }
    for (int i = 0; i < params["numI1"]; i++) {
      _undiagnosed1.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("I1");
    }
    for (int i = 0; i < params["numI2"]; i++) {
      _undiagnosed2.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("I2");
    }
    for (int i = 0; i < params["numJ0"]; i++) {
      _diagnosed0.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("J0");
    }
    for (int i = 0; i < params["numJ1"]; i++) {
      _diagnosed1.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("J1");
    }
    for (int i = 0; i < params["numJ2"]; i++) {
      _diagnosed2.push_back(insert_lineage(params["polytomy_time"]));
      _states.push_back("J2");
    }
    _I0.push_back(params["numI0"]);
    _I1.push_back(params["numI1"]);
    _I2.push_back(params["numI2"]);
    _J0.push_back(params["numJ0"]);
    _J1.push_back(params["numJ1"]);           
    _J2.push_back(params["numJ2"]);           
    _times.push_back(0.0);
    _I0_pop_integral = 0;
    _I1_pop_integral = 0;
    _I2_pop_integral = 0;
  }

  //rMeasure
  string rmeasure(map<string,double> & params, double sample_time, int nlocus)
  {
    throw runtime_error("rmeasure incorporated into simulate function");
  }
	
  //dMeasure
  double dmeasure(map<string,double> & params, 
		  const vector<string> & seqs,
		  int seq_index, 
		  double sample_time) 
  {
    // Partition the likelihood into to two components
    double seq_cll = 0, sample_cll = 0;

    // Extract counts of undiagnosed individuals in each stage
    double numI0 = (double) _undiagnosed0.size();
    double numI1 = (double) _undiagnosed1.size();
    double numI2 = (double) _undiagnosed2.size();

    // If there are no individuals to diagnose, return -INF
    if(numI0 + numI1 + numI2 == 0) {
      return NEG_INF;
    }

    // Randomly choose an undiagnosed individual to diagnose
    // Attach a sequence if there is one
    double total, wI0, wI1, wI2, u;
    wI0  = numI0 * params["diagI0"];
    wI1 =  numI1 * params["diagI1"];
    wI2 =  numI2 * params["diagI2"]; 
    total = wI0 + wI1 + wI2;
    u = runif(0,total);
    int who;

    if(u < wI0) {
      who = (int) floor(runif(0, numI0));
      _diagnosed0.push_back(_undiagnosed0[who]);
      remove_meristem(_undiagnosed0, who);
      if(seqs[seq_index] != "NA") seq_cll = attach_sequence(_diagnosed0.back(),
							   seqs,
							   sample_time,
							   seq_index,
							   params);      
    } else if (u < wI0 + wI1) {
      who = (int) floor(runif(0, numI1));
      _diagnosed1.push_back(_undiagnosed1[who]);
      remove_meristem(_undiagnosed1, who);
      if(seqs[seq_index] != "NA") seq_cll = attach_sequence(_diagnosed1.back(),
							   seqs,
							   sample_time,
							   seq_index,
							   params);      
    } else {
      who = (int) floor(runif(0, numI2));
      _diagnosed2.push_back(_undiagnosed2[who]);
      remove_meristem(_undiagnosed2, who);
      if(seqs[seq_index] != "NA") seq_cll = attach_sequence(_diagnosed2.back(),
							   seqs,
							   sample_time,
							   seq_index,
							   params);      
    }	
    
    // Compute the likelihood of observing the sample
    sample_cll = -(params["diagI0"] * _I0_pop_integral + 
		   params["diagI1"] * _I1_pop_integral + 
		   params["diagI2"] * _I2_pop_integral) + log(total);

    // if(sample_cll > 1000){
      
    //   cout << "large diagnosis likelihood" << endl;
    //   cout << "diagnosis rate I0: " << params["diagI0"] <<  endl;
    //   cout << "diagnosis rate I1: " << params["diagI1"] <<  endl;
    //   cout << "diagnosis rate I2: " << params["diagI2"] <<  endl;
    //   cout << "I0 pop integral: " << _I0_pop_integral <<  endl;      
    //   cout << "I1 pop integral: " << _I1_pop_integral <<  endl;      
    //   cout << "I2 pop integral: " << _I2_pop_integral <<  endl;      
    //   cout << "total: " << total << endl;
    //   cout << "number of I0: " << numI0 << endl;
    //   cout << "number of I1: " << numI1 << endl;
    //   cout << "number of I2: " << numI2 << endl;
    // }
    
    // if(seq_cll > 0){
    //   cout << "large seq probability" << endl;
    //   cout << seq_cll << endl;
    // }

    // Reset the popIntegrals
    _I0_pop_integral = 0;
    _I1_pop_integral = 0;
    _I2_pop_integral = 0;

    //Return the conditional log likelihood
    return sample_cll + seq_cll;
  }


  // Print quantities in the diagnosis likelihood to file
  // void print_particle(const char * particle_file, map<string,double> & params){
  //   // Extract counts of undiagnosed individuals in each stage
  //   double numI0 = (double) _undiagnosed0.size();
  //   double numI1 = (double) _undiagnosed1.size();
  //   double numI2 = (double) _undiagnosed2.size();
    
  //   ofstream filestream;
  //   filestream.open(particle_file);
  //   filestream << "diagnosis rate I0: " << params["diagI0"] <<  endl;
  //   filestream << "diagnosis rate I1: " << params["diagI1"] <<  endl;
  //   filestream << "diagnosis rate I2: " << params["diagI2"] <<  endl;
  //   filestream << "I0 pop integral: " << _I0_pop_integral <<  endl;      
  //   filestream << "I1 pop integral: " << _I1_pop_integral <<  endl;      
  //   filestream << "I2 pop integral: " << _I2_pop_integral <<  endl;      
  //   filestream << "number of I0: " << numI0 << endl;
  //   filestream << "number of I1: " << numI1 << endl;
  //   filestream << "number of I2: " << numI2 << endl;
  //   filestream.close();
  // }


  // Function to compute the log of the hazard of a diagnosis
  // NOTE: to be accurate, this function must be called before dmeasure, which 
  //       modifies the counts of individuals in each class
  double log_diagnosis_hazard(map<string,double> & params){
    // Extract counts of undiagnosed individuals in each stage
    double numI0 = (double) _undiagnosed0.size();
    double numI1 = (double) _undiagnosed1.size();
    double numI2 = (double) _undiagnosed2.size();
    double total, wI0, wI1, wI2;
    wI0 = numI0 * params["diagI0"];
    wI1 = numI1 * params["diagI1"];
    wI2 = numI2 * params["diagI2"]; 
    total = wI0 + wI1 + wI2;
    if(total == 0) return NEG_INF;
    return log(total);
  }

  // Function to compute the probability of no diagnosis
  // NOTE: to be accurate, this function must be called before dmeasure, which 
  //       sets the population integrals to zero
  double log_prob_no_diagnosis(map<string,double> & params){
   return -(params["diagI0"] * _I0_pop_integral + 
	    params["diagI1"] * _I1_pop_integral + 
	    params["diagI2"] * _I2_pop_integral);
  }

  //rProcess
  void rprocess(map<string,double> & params, double start_time, double end_time)
    
  {
    double t = start_time;
    double dt, u; 
    int who;
    double emig = params["emigration_rate"];
    // Counts of individuals in each class
    double numI0, numI1, numI2, numJ0, numJ1, numJ2; 

    while(TRUE)
      {

	// Extract counts of individuals in each class
	numI0 = (double) _undiagnosed0.size();
	numI1 = (double) _undiagnosed1.size();
	numI2 = (double) _undiagnosed2.size();
	numJ0 = (double) _diagnosed0.size();
	numJ1 = (double) _diagnosed1.size();
	numJ2 = (double) _diagnosed2.size();

	// Simulate the time to the next event
	int nrates = 17;
	double rates[nrates];
	rates[0] = params["epsilonI0"]* numI0;
	rates[1] = params["epsilonI1"]* numI1;
	rates[2] = params["epsilonI2"]* numI2;
	rates[3] = params["epsilonJ0"]* numJ0;
	rates[4] = params["epsilonJ1"]* numJ1;
	rates[5] = params["epsilonJ2"]* numJ2;
	rates[6] = params["immigration_rate"];
	rates[7] = (params["muI0"] + emig) * numI0;
	rates[8] = (params["muI1"] + emig) * numI1; 
	rates[9] = (params["muI2"] + emig) * numI2;
	rates[10] = (params["muJ0"] + emig) * numJ0;
	rates[11] = (params["muJ1"] + emig) * numJ1;
	rates[12] = (params["muJ2"] + emig) * numJ2;
	rates[13] = params["gammaI0I1"] * numI0;
	rates[14] = params["gammaI1I2"] * numI1;
	rates[15] = params["gammaJ0J1"] * numJ0;
	rates[16] = params["gammaJ1J2"] * numJ1;
	// Convert to cumulative rates
	for(int i = 1; i < nrates; i++) rates[i] += rates[i-1]; 
	dt = rexp(1/rates[nrates-1]);
	// dt = -log(runif())/rates[nrates-1]; // AAK- CHECK TO SEE IF REXP IS WHAT YOU THINK IT IS.

	// Check the time, leave loop if time has exceeded end_time
	if(t + dt > end_time || !isfinite(dt)) {
	  // Before exiting, add last bit to the popIntegrals
	  _I0_pop_integral += numI0 * (end_time - t);
	  _I1_pop_integral += numI1 * (end_time - t);
	  _I2_pop_integral += numI2 * (end_time - t);
	  break;	
	}

	// Update time
	t += dt;

	// Accumulate the integrals of I0 and I1
	_I0_pop_integral += numI0 * dt;
	_I1_pop_integral += numI1 * dt;
	_I2_pop_integral += numI2 * dt;

	// Determine the type of the event 
	u = runif(0,rates[nrates-1]);
	int eventType = 0;
	while(u > rates[eventType]) eventType++;

	// Change the state of the system according to what event occurred
	switch(eventType) 
	  
	  {
	    
	  case 0:  // Infection event by I0
	    
	    {
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed0[who], t));
	      break;
	    }

	  case 1:  // Infection event by I1
	    
	    {
	      who = (int) floor(runif(0, _undiagnosed1.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed1[who], t));
	      break;
	    }

	  case 2:  // Infection by individual of class I2

	    {
	      who = (int) floor(runif(0, _undiagnosed2.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed2[who], t));
	      break;
	    }	      

	  case 3:  // Infection by individual of class J0
	    
	    {
	      who = (int) floor(runif(0, _diagnosed0.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed0[who], t));
	      break;
	    }

	  case 4:  // Infection by individual of class J1
	    
	    {
	      who = (int) floor(runif(0, _diagnosed1.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed1[who], t));
	      break;
	    }

	  case 5:  // Infection by individual of class J2
	    
	    {
	      who = (int) floor(runif(0, _diagnosed2.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed2[who], t));
	      break;
	    }

	  case 6:  // Immigration event 
	    
	    {
	      // Note that here we are assuming that immigration means being infected with 
	      // a distantly related sequence
	      _undiagnosed0.push_back(insert_lineage(params["polytomy_time"]));
	      break;
	    }

	  case 7: // Death of an I0 individual

	    {	    
	      // Randomly choose an I0 victim
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      terminate_meristem(_undiagnosed0[who], t);						
	      remove_meristem(_undiagnosed0, who);
	      break;
	    }
	    
	  case 8: // Death of an I1 individual
	    
	    {
	      // Randomly choose an I1 victim
	      who = (int) floor(runif(0, _undiagnosed1.size())); 
	      terminate_meristem(_undiagnosed1[who], t);						
	      remove_meristem(_undiagnosed1, who);
	      break;
	    }

	  case 9: // Death of an I2 individual
	    
	    {
	      // Randomly choose an I1 victim
	      who = (int) floor(runif(0, _undiagnosed2.size())); 
	      terminate_meristem(_undiagnosed2[who], t);						
	      remove_meristem(_undiagnosed2, who);
	      break;
	    }

	  case 10: // Death of J0 individual
	    
	    {
	      // Randomly choose a J0 victim
	      who = (int) floor(runif(0,_diagnosed0.size())); 
	      terminate_meristem(_diagnosed0[who],t);						
	      remove_meristem(_diagnosed0,who);
	      break;
	    }

	  case 11: // Death of a J1 individual
	    
	    {
	      // Randomly choose a J1 victim
	      who = (int) floor(runif(0, _diagnosed1.size())); 
	      terminate_meristem(_diagnosed1[who], t);						
	      remove_meristem(_diagnosed1, who);
	      break;
	    }
	    
	  case 12: // Death of a J2 individual
	    
	    {
	      // Randomly choose a J2 victim
	      who = (int) floor(runif(0, _diagnosed2.size())); 
	      terminate_meristem(_diagnosed2[who], t);						
	      remove_meristem(_diagnosed2, who);
	      break;
	    }

	  case 13: // I0 -> I1
	    
	    {
	      // Randomly choose an I0 to transition
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      // Update meristem lists
	      _undiagnosed1.push_back(_undiagnosed0[who]);
	      remove_meristem(_undiagnosed0, who);						
	      break;
	    }

	  case 14: // I1 -> I2
	    
	    {
	      // Randomly choose an I1 to transition
	      who = (int) floor(runif(0, _undiagnosed1.size())); 
	      // Update meristem lists
	      _undiagnosed2.push_back(_undiagnosed1[who]);
	      remove_meristem(_undiagnosed1, who);						
	      break;
	    }	    

	  case 15: // J0 -> J1
	    
	    {
	      // Randomly choose a J0 to transition
	      who = (int) floor(runif(0, _diagnosed0.size())); 
	      // Update meristem lists
	      _diagnosed1.push_back(_diagnosed0[who]);
	      remove_meristem(_diagnosed0,who);			
	      break;
	    }
	    
	  case 16: // J1 -> J2

	    {
	      // Randomly choose a J1 to transition
	      who = (int) floor(runif(0, _diagnosed1.size())); 
	      // Update meristem lists
	      _diagnosed2.push_back(_diagnosed1[who]);
	      remove_meristem(_diagnosed1, who);			
	      break;
	    }
	    
	  default:

	    {
	      throw runtime_error("Impossible event in Gillespie Algorithm!");
	      break;
	    }
	  } 
      }
  }
   		
  //simulate
  void simulate(map<string,double> & params, double start_time, double end_time, string root_sequence)
  {
    
    double t = start_time;
    double dt, u; 
    int who;
    double emig = params["emigration_rate"];
    // Counts of individuals in each class
    double numI0, numI1, numI2, numJ0, numJ1, numJ2; 

    // Size of sequences to simulate
    int nlocus = root_sequence.size();
    // Seed the sequence at the root
    set_root_sequence(root_sequence);

    while(TRUE)
      {
	// Extract counts of individuals in each class
	numI0 = (double) _undiagnosed0.size();
	numI1 = (double) _undiagnosed1.size();
	numI2 = (double) _undiagnosed2.size();
	numJ0 = (double) _diagnosed0.size();
	numJ1 = (double) _diagnosed1.size();
	numJ2 = (double) _diagnosed2.size();
	
	// Simulate the time to the next event
	int nrates = 20;
	double rates[nrates];
	rates[0] = params["epsilonI0"]* numI0;
	rates[1] = params["epsilonI1"]* numI1;
	rates[2] = params["epsilonI2"]* numI2;
	rates[3] = params["epsilonJ0"]* numJ0;
	rates[4] = params["epsilonJ1"]* numJ1;
	rates[5] = params["epsilonJ2"]* numJ2;
	rates[6] = params["immigration_rate"];
	rates[7] = (params["muI0"] + emig) * numI0;
	rates[8] = (params["muI1"] + emig) * numI1; 
	rates[9] = (params["muI2"] + emig) * numI2;
	rates[10] = (params["muJ0"] + emig) * numJ0;
	rates[11] = (params["muJ1"] + emig) * numJ1;
	rates[12] = (params["muJ2"] + emig) * numJ2;
	rates[13] = params["gammaI0I1"] * numI0;
	rates[14] = params["gammaI1I2"] * numI1;
	rates[15] = params["gammaJ0J1"] * numJ0;
	rates[16] = params["gammaJ1J2"] * numJ1;
	rates[17] = params["diagI0"] * numI0;
	rates[18] = params["diagI1"] * numI1;
	rates[19] = params["diagI2"] * numI2;
	// Convert to cumulative rates
	for(int i = 1; i < nrates; i++) rates[i] += rates[i-1]; 
	dt = rexp(1/rates[nrates-1]);
	
	// Update time
	t += dt;
	
	// Check the time, leave loop if time has exceeded end_time
	if(t > end_time || !isfinite(dt)) {
	  // Before leaving the loop, terminate all meristems and 
	  // add them to the tree
	  for (int i = 0; i < numI0; i++) {
	    terminate_meristem(_undiagnosed0.at(i), end_time);
	    _states.push_back("I0");
	  }
	  for (int i = 0; i < numI1; i++) {
	    terminate_meristem(_undiagnosed1.at(i), end_time);
	    _states.push_back("I1");
	  }
	  for (int i = 0; i < numI2; i++) {
	    terminate_meristem(_undiagnosed2.at(i), end_time);
	    _states.push_back("I2");
	  }
	  for (int i = 0; i < numJ0; i++) {
	    terminate_meristem(_diagnosed0.at(i), end_time);
	    _states.push_back("J0");
	  }
	  for (int i = 0; i < numJ1; i++) {
	    terminate_meristem(_diagnosed1.at(i), end_time);
	    _states.push_back("J1");
	  }
	  for (int i = 0; i < numJ2; i++) {
	    terminate_meristem(_diagnosed2.at(i), end_time);
	    _states.push_back("J2");
	  }
	  break;	
	}
	
	// Determine the type of the event 
	u = runif(0,rates[nrates-1]);
	int event_type = 0;
	while(u > rates[event_type]) event_type++;
	
	
	// Change the state of the system according to what event occurred
	switch(event_type) 

	  {

	  case 0 :  // Infection event by I0
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed0[who], t));
	      _states.push_back("I0");
	      // Update counts of individuals
	       update_state_counts(1,0,0,0,0,0);
	       break;
	    }
	    
	  case 1 :  //Infection event by I1
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0, _undiagnosed1.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed1[who], t));
	      _states.push_back("I1");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }
	    
	  case 2 :  //Infection by individual of class I2
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0,_undiagnosed2.size())); 
	      _undiagnosed0.push_back(add_meristem(_undiagnosed2[who], t));
	      _states.push_back("I2");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }	      
	    
	  case 3 :  //Infection by individual of class J0
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0,_diagnosed0.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed0[who], t));
	      _states.push_back("J0");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }
	    
	  case 4 :  //Infection by individual of class J1
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0,_diagnosed1.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed1[who], t));
	      _states.push_back("J1");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }
	    
	  case 5 :  //Infection by individual of class J2
	    
	    {
	      // Choose random individual as infector and modify tree
	      who = (int) floor(runif(0,_diagnosed2.size())); 
	      _undiagnosed0.push_back(add_meristem(_diagnosed2[who], t));
	      _states.push_back("J2");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }
	    
	  case 6 :  //Immigration event 
	    
	    {
	      // Note that here we are assuming that immigration means being infected with 
	      // a distantly related sequence (and so a new individual enters the I0 class
	      _undiagnosed0.push_back(insert_lineage(params["polytomy_time"]));
	      _states.push_back("I0");
	      // Update counts of individuals
	      update_state_counts(1,0,0,0,0,0);
	      break;
	    }
	    
	  case 7 : // Death of an I0 individual
	    
	    {	    
	      // Randomly choose an I0 victim
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      terminate_meristem(_undiagnosed0[who], t);						
	      remove_meristem(_undiagnosed0, who);
	      _states.push_back("I0");	     
	      // Update counts of individuals
	      update_state_counts(-1,0,0,0,0,0);
	      break;
	    }
	    
	  case 8 : //Death of an I1 individual
	    
	    {
	      // Randomly choose an I1 victim
	      who = (int) floor(runif(0, _undiagnosed1.size())); 
	      terminate_meristem(_undiagnosed1[who], t);						
	      remove_meristem(_undiagnosed1, who);
	      _states.push_back("I1");	     
	      // Update counts of individuals
	      update_state_counts(0,-1,0,0,0,0);
	      break;
	    }
	    
	  case 9 : //Death of an I2 individual
	    
	    {
	      // Randomly choose an I1 victim
	      who = (int) floor(runif(0, _undiagnosed2.size())); 
	      terminate_meristem(_undiagnosed2[who], t);						
	      remove_meristem(_undiagnosed2, who);
	      _states.push_back("I2");	     
	      // Update counts of individuals
	      update_state_counts(0,0,-1,0,0,0);
	      break;
	    }
	    
	  case 10 : //Death of J0 individual
	    
	    {
	      // Randomly choose a J0 victim
	      who = (int) floor(runif(0, _diagnosed0.size())); 
	      terminate_meristem(_diagnosed0[who], t);						
	      remove_meristem(_diagnosed0, who);
	      _states.push_back("J0");	     
	      // Update counts of individuals
	      update_state_counts(0,0,0,-1,0,0);
	      break;
	    }
	    
	  case 11 : //Death of a J1 individual
	    
	    {
	      // Randomly choose a J1 victim
	      who = (int) floor(runif(0, _diagnosed1.size())); 
	      terminate_meristem(_diagnosed1[who], t);						
	      remove_meristem(_diagnosed1, who);
	      _states.push_back("J1");	     
	      // Update counts of individuals
	      update_state_counts(0,0,0,0,-1,0);
	      break;
	    }
	    
	  case 12 : //Death of a J2 individual
	    
	    {
	      // Randomly choose a J1 victim
	      who = (int) floor(runif(0, _diagnosed2.size())); 
	      terminate_meristem(_diagnosed2[who], t);						
	      remove_meristem(_diagnosed2, who);
	      _states.push_back("J2");	     
	      // Update counts of individuals
	      update_state_counts(0,0,0,0,0,-1);
	      break;
	    }
	    
	  case 13 : //I0 -> I1
	    
	    {
	      // Randomly choose an I0 to transition
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      //Branch the lineage so that we can track the state transition (but don't note the new branch)
	      int ignore = add_leaf(_undiagnosed0[who], t, t);	   
	      _states.push_back("I0");	   
	      // Update meristem lists
	      _undiagnosed1.push_back(_undiagnosed0[who]);
	      remove_meristem(_undiagnosed0, who);						
	      // Update counts of individuals
	      update_state_counts(-1,1,0,0,0,0);
	      break;
	    }
	    
	  case 14 : //I1 -> I2
	    
	    {
	      // Randomly choose an I1 to transition
	      who = (int) floor(runif(0,_undiagnosed1.size())); 
	      //Branch the lineage so that we can track the state transition (but don't note the new branch)
	      int ignore = add_leaf(_undiagnosed1[who], t, t);	   
	      _states.push_back("I1");	   
	      // Update meristem lists
	      _undiagnosed2.push_back(_undiagnosed1[who]);
	      remove_meristem(_undiagnosed1, who);						
	      // Update counts of individuals
	      update_state_counts(0,-1,1,0,0,0);
	      break;
	    }	    

	  case 15 : //J0 -> J1
	    
	    {
	      // Randomly choose a J0 to transition
	      who = (int) floor(runif(0,_diagnosed0.size())); 
	      //Branch the lineage so that we can track the state transition (but don't note the new branch)
	      int ignore = add_leaf(_diagnosed0[who], t, t);	   
	      _states.push_back("J0");	   
	      // Update meristem lists
	      _diagnosed1.push_back(_diagnosed0[who]);
	      remove_meristem(_diagnosed0, who);			
	      // Update counts of individuals
	      update_state_counts(0,0,0,-1,1,0);
	      break;
	    }
	    
	  case 16 : // J1 -> J2

	    {
	      // Randomly choose a J1 to transition
	      who = (int) floor(runif(0,_diagnosed1.size())); 
	      //Branch the lineage so that we can track the state transition (but don't note the new branch)
	      int ignore = add_leaf(_diagnosed1[who], t, t);	   
	      _states.push_back("J1");	   
	      // Update meristem lists
	      _diagnosed2.push_back(_diagnosed1[who]);
	      remove_meristem(_diagnosed1, who);			
	      // Update counts of individuals
	      update_state_counts(0,0,0,0,-1,1);
	      break;
	    }

	     
	  case 17: // I0 -> J0

	    {
	      // Randomly choose an I0 to diagnose
	      who = (int) floor(runif(0, _undiagnosed0.size())); 
	      // Update meristem lists
	      _diagnosed0.push_back(_undiagnosed0[who]);
	      remove_meristem(_undiagnosed0, who);						
	      //Simulate a sample for that diagnosed individual with probability psequence
	      double p = runif(0,1);
	      if(p < params["psequence"]) {
		_seqs.push_back(simulate_sequence(_diagnosed0.back(),
						  nlocus,
						  t,
						  params));
		_sequence_times.push_back(t);
		_states.push_back("I0"); // one for the call to add_meristem	   
		_states.push_back("Sequence"); // one for the call to terminate_meristem	   	     
	      } else {
		//Diagnosis only 
		_seqs.push_back("NA");
		_sequence_times.push_back(t);
		int ignore = add_meristem(_diagnosed0.back(), t);
		_states.push_back("I0"); // one for the call to add_meristem only	   	     
	      }
	      // Update counts of individuals
	      update_state_counts(-1,0,0,1,0,0);
	      break;
	    }
	     	     
	    
	  case 18: //I1 -> J1
	    
	    {
	      // Randomly choose an I1 to diagnose
	      who = (int) floor(runif(0,_undiagnosed1.size())); 
	      // Update meristem lists
	      _diagnosed1.push_back(_undiagnosed1[who]);
	      remove_meristem(_undiagnosed1, who);						
	      //Simulate a sample for that diagnosed individual with probability psequence
	      double p = runif(0,1);
	      if(p < params["psequence"]) {
		_seqs.push_back(simulate_sequence(_diagnosed1.back(),
						  nlocus,
						  t,
						  params));
		_sequence_times.push_back(t);
		_states.push_back("I1"); // one for the call to add_meristem	   
		_states.push_back("Sequence"); // one for the call to terminate_meristem	   	     
	      } else {
		//Diagnosis only 
		_seqs.push_back("NA");
		_sequence_times.push_back(t);
		int ignore = add_meristem(_diagnosed1.back(), t);
		_states.push_back("I1"); // one for the call to add_meristem only	   	     
	      }
	      // Update counts of individuals
	      update_state_counts(0,-1,0,0,1,0);
	      break;
	    }
	    
	  case 19: //I2 -> J2
	    
	    {
	      // Randomly choose an I2 to diagnose
	      who = (int) floor(runif(0, _undiagnosed2.size())); 
	      // Update meristem lists
	      _diagnosed2.push_back(_undiagnosed2[who]);
	      remove_meristem(_undiagnosed2, who);						
	      //Simulate a sample for that diagnosed individual with probability psequence
	      double p = runif(0,1);
	      if(p < params["psequence"]) {
		_seqs.push_back(simulate_sequence(_diagnosed2.back(),
						  nlocus,
						  t,
						  params));
		_sequence_times.push_back(t);
		_states.push_back("I2"); // one for the call to add_meristem	   
		_states.push_back("Sequence"); // one for the call to terminate_meristem	   	     
	      } else {
		// Diagnosis only 
		_seqs.push_back("NA");
		_sequence_times.push_back(t);
		int ignore = add_meristem(_diagnosed2.back(), t);
		_states.push_back("I2"); // one for the call to add_meristem only	   	     
	      }
	      // Update counts of individuals
	      update_state_counts(0,0,-1,0,0,1);
	      break;
	    }
	    
	  default :

	    {
	      throw runtime_error("Impossible event in Gillespie Algorithm!");
	      break;
	    }
	  } 
	_times.push_back(t);	   
      }
  }
  
  //Function to remove a meristem from the vector of meristems
  //This function was written as is to avoid reallocation of the vector when
  // deleting a single element
  void remove_meristem (vector<node_index> & m, int victim_index)
  {
    //Copy last element in the meristem vector over the meristem at victim_index
    m.at(victim_index) = m[m.size() - 1];
    //Delete last element in the meristem vector
    m.erase(--m.end());	
  }	

  //Function to append updated counts of individuals in each state onto the state vectors
  void update_state_counts (int dI0, int dI1, int dI2, int dJ0, int dJ1, int dJ2) {
	  _I0.push_back(_I0.back() + dI0);
	  _I1.push_back(_I1.back() + dI1);
	  _I2.push_back(_I2.back() + dI2);	   
	  _J0.push_back(_J0.back() + dJ0);
	  _J1.push_back(_J1.back() + dJ1);          
	  _J2.push_back(_J2.back() + dJ2);                
  }
  
  //Function to print sequences
  void printSeqs() {
    for(int i = 0; i < _seqs.size(); i++)
      cout << _sequence_times.at(i) << ' ' << _seqs.at(i) << endl;
  }

  //Function to save sequences to a file
  void save_seqs(char *seqfile) {
    ofstream fileStream;
    fileStream.open(seqfile);
    //fileStream << "time seq" << endl;
    for(int i = 0; i < _seqs.size(); i++)
      fileStream << _sequence_times.at(i) << ' ' << _seqs.at(i) << endl;
    fileStream.close();
  }

  //Function to print counts of individuals through time
  void printCounts() {
    cout << "time I0 I1 J0 J1" << endl;
    for(int i = 0; i < _times.size(); i++)
      cout << _times.at(i) << ' ' 
	   << _I0.at(i) << ' ' 
	   << _I1.at(i) << ' ' 
	   << _J0.at(i) << ' ' 
	   << _J1.at(i) << endl;
  }

//80////////////////////////////////////////////////////////////////////////////

  //Function to save counts of individuals through time to a file
  void save_counts(char *popfile) {
    ofstream fileStream;
    fileStream.open(popfile);
    fileStream << "time I0 I1 I2 J0 J1 J2" << endl;
    for(int i = 0; i < _times.size(); i++)
      fileStream << _times.at(i) << ' ' << _I0.at(i) << ' ' 
		 << _I1.at(i) << ' ' <<  _I2.at(i) << ' ' 
		 << _J0.at(i) << ' ' << _J1.at(i) << ' ' 
		 <<  _J2.at(i) << endl;
    fileStream.close();
  }

  // Sets up the header in the file that stores the state counts
  void write_states_header(const char *states_file){
    ofstream filestream;
    filestream.open(states_file);
    filestream << "data_point particle_index I0 I1 I2 J0 J1 J2" << endl;
    filestream.close();
  }

  // Function to save counts of individuals at a specific data point
  // Intended to be called while running a particle filter
  // NOTE: the function 'write_states_header' should be called before this one
  void write_states(int data_point, int particle_index, const char *states_file){
    ofstream file_stream;
    file_stream.open(states_file, std::ios::app);
    file_stream << data_point << ' ' << particle_index << ' ';
    file_stream << _undiagnosed0.size() << ' ' << _undiagnosed1.size() << ' ';
    file_stream << _undiagnosed2.size() << ' ' << _diagnosed0.size() << ' ';
    file_stream << _diagnosed1.size() << ' ' << _diagnosed2.size() << endl;
    file_stream.close();
  }

  //Function to print ouch tree
  void printOuchTree() {
    cout << "node ancestor time state" << endl;
    for(int i = 0; i < _nodes.size(); i ++) {
      if(i == 0) {
	cout << _nodes.at(i) << " NA " << _node_times.at(i) << ' ' 
	     << _states.at(i) << endl;   
      } else {
	cout << _nodes.at(i) << ' ' << _ancestors.at(i) << ' ' 
	     << _node_times.at(i) << ' ' << _states.at(i) << endl;   
      }
    }
  }
  
  //Function to save tree table to a file
  void save_tree_table(char *tree_file) 
  {
    ofstream fileStream;
    fileStream.open(tree_file);
    fileStream << "node ancestor id time state" << endl;
    for(int i = 0; i < _nodes.size(); i ++) {
      if(i == 0) {
	fileStream << _nodes.at(i) << " NA " 
		   << _ids.at(i) << ' ' 
		   << _node_times.at(i) << ' ' 
		   << _states.at(i) << endl;   
      } else {
	fileStream << _nodes.at(i) << ' ' 
		   << _ancestors.at(i) << ' ' 
		   << _ids.at(i) << ' ' 
		   << _node_times.at(i) << ' ' 
		   << _states.at(i) << endl;   
      }
    }
    fileStream.close();
  }
  
};

//80////////////////////////////////////////////////////////////////////////////

#endif

