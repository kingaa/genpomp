// -*- mode: C++; -*-
// This is the implementation file for the class basenode

#include "basenode.h"

//Default Constructor
basenode::basenode (void) : _mother(MOTHER_EVE), _time_stamp(NAN) { }

//Constructor
basenode::basenode (node_index mom, double t) : 
  _mother(mom), _time_stamp(t) { }

//Get functions
node_index basenode::mother (void) const {
  return _mother;
}
double basenode::time (void) const {
  return _time_stamp;
}

//Set functions
void basenode::mother (node_index mom) {
  _mother = mom;
}
void basenode::time (double t) {
  _time_stamp = t;
}

//Check functions
bool basenode::is_meristem (void) const {
  return isnan(_time_stamp);
}
