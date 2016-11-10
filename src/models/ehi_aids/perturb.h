// -*- mode: C++; -*-
#ifndef _PERTURB_H_
#define _PERTURB_H_

#include "abstractPerturb.h"
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <Rmath.h>

using namespace std;
typedef map<string, double>::iterator iter_type;

class Perturb : public AbstractPerturb
{

 private :

  // Random walk standard deviations for regular parameters
  map<string, double> rw_sds;
  // Random walk standard deviations for IVPs
  map<string, double> ivp_rw_sds; 

 public :
  
  // Constructor
  Perturb(map<string, double> random_walk_sds, 
	  map<string, double> ivp_random_walk_sds) {
    rw_sds = random_walk_sds;
    ivp_rw_sds = ivp_random_walk_sds;
  }

  // Function to perturb parameters (written by the user)
  void perturb_parameters(map<string, double> & params, double alpha) {
    // Perturb all 'independent' parameters (that is, parameters that do not
    // need to be perturbed jointly
    double sd;
    for(iter_type iterator = rw_sds.begin(); iterator != rw_sds.end(); iterator++) {    
      // Cool the standard deviation
      sd = alpha * iterator->second;
      // Perturb and write over old parameter value
      params[iterator->first] *= rlnorm(-(sd * sd)/2,sd);
    }
  }
  
  //Function to perturb inital value parameters (written by the user)
  void perturb_initial_value_parameters(map<string, double> & params, double alpha){
    // Perturb the initial value parameters
    double sd;
    for(iter_type iterator = ivp_rw_sds.begin(); iterator != ivp_rw_sds.end(); iterator++) {    
      // Cool the standard deviation
      sd = alpha * iterator->second;
      // Perturb and write over old parameter value
      params[iterator->first] *= rlnorm(-(sd * sd)/2,sd);
    }
  }

};

#endif
