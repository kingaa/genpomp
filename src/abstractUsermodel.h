// -*- C++ -*-
//Abstract base class to prescribe interface for usermodels

#ifndef _ABSTRACT_USERMODEL_H_
#define _ABSTRACT_USERMODEL_H_

#include "type_defs.h"
#include "gsl_rng.h"

using namespace std;

class AbstractUsermodel : public tree
{

public:

  //constructor
  AbstractUsermodel(substModel model, 
		    map<string,double> & params, 
		    bool make_tree_table,
		    int num_gnodes = 0) : tree(model,
					       params["polytomy_time"],
					       make_tree_table,
					       num_gnodes) {}
  
  //rMeasure
  virtual string rmeasure(map<string,double> & params, 
			  double sample_time, 
			  int nlocus) = 0;
	
  //dMeasure
  virtual double dmeasure(gsl_rng * rngptr,
			  map<string,double> & params, 
			  const vector<string> & seqs,
			  int seq_index, 
			  double sample_time) = 0;
		
  //rProcess
  virtual void rprocess(gsl_rng * rngptr,
			map<string,double> & params, 
			double start_time, 
			double end_time) = 0;
};

#endif
