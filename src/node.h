// -*- mode: C++; -*-
#ifndef _NODE_H_
#define _NODE_H_

#include <cmath>
#include "basenode.h"
#include "gnode.h"
#include "type_defs.h"

#define NO_LINEAGE (-1)

using namespace std;

class node : public basenode 
{
private:
  gnode_index _lineage; // a breadcrumb that indicates the genetic lineage that this node belongs to

public:
  node (void);
  node (node_index mom, double t = NAN, gnode_index lin = NO_LINEAGE);
  gnode_index get_lineage(void) const;
  void set_lineage (gnode_index gptr);
};

#endif
