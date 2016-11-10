// -*- mode: C++; -*-/
//Header file for basenode class
//This is primarily intended as an ancestor class for node types to
//	represent the infection and genetic trees

#ifndef _BASENODE_H_
#define _BASENODE_H_

#include <cmath>
#include "type_defs.h"

#define MOTHER_EVE (-1)

class basenode {

private:
  node_index _mother;
  double _time_stamp;

public:
  basenode (void);	//Default Constructor
  basenode (node_index mom, double current_time = NAN);	//Constructor
  node_index mother (void) const; //Get functions
  double time (void) const;
  void mother (node_index mom); //Set functions
  void time (double t);
  bool is_meristem (void) const; //Check functions
  
};

#endif
