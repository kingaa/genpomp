// -*- mode: C++; -*-
/*
Node.cpp
Implementation file for the node class
*/

#include "node.h"

node::node (void) : basenode(), _lineage(NO_LINEAGE) { }

node::node (node_index mom, double t, gnode_index lin) : 
	basenode(mom,t), _lineage(lin) { }

gnode_index node::get_lineage(void) const{
	return _lineage;
}

void node::set_lineage (gnode_index new_lineage) {
	_lineage = new_lineage;
}
