// -*- mode: C++; -*-
/*
gnode.cpp
Implementations of member functions of the gnode class
*/

#include "gnode.h"

gnode::gnode (void) : basenode(), _daughter1(NO_DAUGHTER),
		      _daughter2(NO_DAUGHTER), _ellmat(Lmatrix()), _evo_branch_length(-1){ }
	
gnode::gnode (node_index mom, double t, node_index d1, 
	      node_index d2, double evo_branch_length): basenode(mom,t), _daughter1(d1),
		     _daughter2(d2), _ellmat(Lmatrix()), _evo_branch_length(evo_branch_length) { }

node_index gnode::daughter_one (void) const {
	return _daughter1;
}

node_index gnode::daughter_two (void) const{
	return _daughter2;
}

void gnode::daughter_one (node_index d1){
	_daughter1 = d1;
}

void gnode::daughter_two (node_index d2){
	_daughter2 = d2;
}

Lmatrix *gnode::ellmatrix (void) {
  return &_ellmat;
}

void gnode::ellmatrix (Lmatrix ell) {
  _ellmat = ell;
}

bool gnode::has_data (void) const {
    return (_ellmat.data != NULL);
  }

void gnode::evo_branch_length(double evo_branch_length){
  _evo_branch_length = evo_branch_length;
}

double gnode::evo_branch_length(){
  return _evo_branch_length;
}
