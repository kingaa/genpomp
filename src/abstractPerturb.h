// -*- mode: C++; -*-
#ifndef _ABSTRACT_PERTURB_H_
#define _ABSTRACT_PERTURB_H_

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <Rmath.h>
#include "type_defs.h"

using namespace std;

class AbstractPerturb {

public :

  // Constructor
  //AbstractPerturb(map<string, double> rw_params) = 0;

  //Function to perturb parameters to be written by the user
  virtual void perturb_parameters(map<string, double> & params, double alpha) = 0;

};

#endif
