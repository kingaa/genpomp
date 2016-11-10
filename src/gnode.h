// -*- mode: C++; -*-
#ifndef _GNODE_H_
#define _GNODE_H_

#include "basenode.h"
#include "lmatrix.h"
#include "substmodel.h"
#include "type_defs.h"

using namespace std;

#define NO_DAUGHTER (-1)

class gnode : public basenode 
{
private:
  node_index _daughter1;
  node_index _daughter2;
  Lmatrix _ellmat; 
  double _evo_branch_length;

public:
  gnode (void);
  gnode (node_index mom, 
	 double t, 
	 node_index d1 = NO_DAUGHTER, 
	 node_index d2 = NO_DAUGHTER, 
	 double evo_branch_length = -1);
  node_index daughter_one (void) const;
  node_index daughter_two (void) const;
  void daughter_one (node_index d1); 
  void daughter_two (node_index d2);   
  Lmatrix *ellmatrix (void);
  void ellmatrix (Lmatrix ell);
  bool has_data (void) const;
  void evo_branch_length(double evo_branch_length);
  double evo_branch_length();
};

#endif
