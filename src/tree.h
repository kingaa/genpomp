// -*- mode: C++; -*-
#ifndef _TREE_H_
#define _TREE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <stdexcept>
#include <omp.h>
#include "lmatrix.h"
#include "basenode.h"
#include "node.h"
#include "gnode.h"
#include "substmodel.h"
#include "type_defs.h"
#include "gsl_rng.h"
#include "gsl_randist.h"

#define NO_SEQ (-1)

using namespace std;

//typedef std::vector<double> doubvec;

class tree : public vector<node>
{

private:	

  vector<gnode> _gtree;
  node_index _dummy_root;
  substModel _model;
  double _depth;
  double _loglik;
  bool _make_tree_table;  // Whether or not to make the tree table 

public:
  
  // Columns of the tree table
  vector<int> _nodes;        
  vector<int> _ancestors;
  vector<double> _node_times;
  vector<int> _ids;

  // Default constructor
  tree (void) : _gtree(), 
		_dummy_root(),
		_model(), 
		_depth(0), 
		_loglik(0),
		_make_tree_table(FALSE),
		_nodes(),
		_ancestors(),
		_node_times(),
		_ids() {
    // Set the first node in the gene tree
    _gtree.push_back(gnode(MOTHER_EVE, 0));
  }

  // Constructor
  tree (substModel model, 
	double polytomy_time,
	bool make_tree_table = FALSE,
	int num_gnodes = 0)  : _gtree(num_gnodes),
			       _dummy_root(),
			       _model(model),
			       _depth(0), 
			       _loglik(0), 
			       _make_tree_table(make_tree_table),
			       _nodes(),
			       _ancestors(),
			       _node_times(),
			       _ids() {
    // Set the first node in the gene tree
    _gtree.push_back(gnode(MOTHER_EVE, polytomy_time));
  }
					
  void set_mol_model(substModel & newModel) {
    _model = newModel;
  }

  double loglik (void) const {
    return _loglik;
  }

  bool is_root (node_index & n) {
    return (at(n).mother() == MOTHER_EVE);
  }
  
  bool is_groot (node_index & n) { 
    return (_gtree.at(n).mother() == MOTHER_EVE);
  }  
  
  bool is_groot (vector<gnode> & gtree, node_index & n) {
    return (gtree.at(n).mother() == MOTHER_EVE);
  }  

  /*
    INPUT: (1) Last fixed node added
           (2) The ancestor of that fixed node
	   (3) The index of the meristem that this node belongs to
	   (4) The time of the fixed node
    MODIFIES: Pushes these values to vectors that will become columns of the 
              tree table.
    OUTPUT: none
  */
  void add_tree_table_row(int node, 
			  int ancestor, 
			  int id,
			  double node_time)
  {
    _nodes.push_back(node);
    _ancestors.push_back(ancestor);
    _ids.push_back(id);
    _node_times.push_back(node_time);
  }

  /*
    INPUT: Time of the polytomy
    MODIFIES: Adds a new meristem, and a new fixed node at the polytomy. 
              Updates _dummy_root to be the index of this new fixed node. This
	      is where the next founder or immigrant will be attached if there
	      is one. 
	      Also updates the vectors that make up the tree table if requested.
    OUTPUT: The index in the vector of meristems of the new meristem.
  */
  // Adds a new meristem with a root at the polytomy
  meristem_index insert_lineage(double polytomy_time) {

    if(size() == 0){

      //First root to be added
      push_back(node(MOTHER_EVE, polytomy_time, NO_LINEAGE)); // true root
      push_back(node(0, polytomy_time, NO_LINEAGE)); //dummy root
      _dummy_root = 1;
      push_back(node(0, polytomy_time, NO_LINEAGE)); // meristem

      // Push states if requested
      if(_make_tree_table) 
	add_tree_table_row(0, MOTHER_EVE, 2, polytomy_time);

      return 2;

    } else {

      // Otherwise attach the new lineage to the current dummy root
      // And attach a new dummy root to the current dummy root
      push_back(node(_dummy_root, polytomy_time, NO_LINEAGE)); // meristem
      node_index new_meristem = size() - 1;
      push_back(node(_dummy_root, polytomy_time, NO_LINEAGE)); // new dummy root
      // Push node to tree table if requested
      if(_make_tree_table) 
	add_tree_table_row(_dummy_root, 
			   at(_dummy_root).mother(), 
			   new_meristem,
			   polytomy_time);

      // Update where to attach new roots
      _dummy_root = size() - 1;
      
      return new_meristem;
    }
  }

  /*
    INPUT: (1) The meristem to branch
           (2) The time of the branch
    MODIFIES: Adds one fixed node to the infection tree at the time of 
              the branch and rewires the meristem to this new mother.
    OUTPUT: Returns the fixed node that is new mother of the meristem.
            The branching process is completed by either add_meristem
	    or add_sequence.
  */
  node_index rewire_meristem (meristem_index me, 
			      double current_time) {
    
    // Get node indices of new mother and new grandma after rewiring
    node_index new_mom = size(); 
    node_index new_grandma = at(me).mother();

    // Rewire, preserving meristem index 
    push_back(node(new_grandma, current_time, NO_LINEAGE)); // add new_mom
    at(me).mother(new_mom);

    // Update depth of tree
    _depth = current_time;

    // Push node to tree table if requested
    if(_make_tree_table) 
      add_tree_table_row(new_mom, new_grandma, me, current_time);

    // Return node_index of new internal fixed node
    return(new_mom);
  }

  /*
    INPUT: (1) The index of the meristem to branch 
           (2) The time to make the branch
    MODIFIES: Adds two nodes: 
      (1) a new fixed node to represent an infection event
      (2) a new meristem to represent a new case
    OUTPUT: index of the new meristem 
  */

  meristem_index add_meristem (meristem_index infector, double current_time) {

    // Rewire meristem to accomodate new branch    
    node_index new_mom = rewire_meristem(infector, current_time);
    // Add new meristem
    push_back(node(new_mom, current_time, NO_LINEAGE));
    //Return index of new meristem
    return size() - 1;
  }

  /*
    INPUT: (1) index of meristem to branch 
           (2) time to make the branch
           (3) time to fix the leaf
    MODIFIES: Adds two fixed nodes to the infection tree: 
      (1) a new internal node
      (2) a new leaf
    OUTPUT: index of the new leaf that is sister to brancher
  */
  node_index add_leaf(meristem_index brancher, 
		      double branch_time,
		      double leaf_time) 
  { 
    // Rewire meristem to accomodate new branch
    node_index new_mom = rewire_meristem(brancher, branch_time);
    // Add new leaf
    push_back(node(new_mom, leaf_time, NO_LINEAGE));
    // Update tree depth
    _depth = branch_time;
    //Return index of new fixed node
    return size() - 1;
  }

  /*
    INPUT: The meristem index to terminate
    MODIFIES: Changes the time of the meristem at new_leaf to 
              current_time. The meaning of this time stamp changes
	      from when to node was created to when it is in time.
	      Also updates the tree table if requested. 
    OUTPUT: none
  */
  void terminate_meristem (meristem_index new_leaf, double current_time) 
  {
    // Fix the time at the meristem
    at(new_leaf).time(current_time);

    // Update tree depth
    _depth = current_time;
    
    // Push node to tree table if requested
    if(_make_tree_table) 
      add_tree_table_row(new_leaf, 
			 at(new_leaf).mother(), 
			 new_leaf,
			 current_time);
  }

  /*
    INPUT: (1) Index where the sequence was attached
           (2) Index of the genetic node associated with this sequence
	   (3) A double to hold the time of mother of this leaf in the genetic
	       tree 
    MODIFIES: Adds breadcrumbs to all nodes in the infection tree that are on 
              the branch subtending the new sequence node. Writes tthe time of 
	      the branch point to (3). Finally, writes over breadcrumbs in the 
	      branching node and those in the nodes on the branch that subtends 
	      that node (if we have not recursed to the root).
    OUTPUT: Returns either:
              (1) The breadcrumb at the branching node that was written over
	      (2) -1, in the case of recursing to the root
  */

  gnode_index update_breadcrumbs(node_index index, 
				 gnode_index new_leaf, 
				 double & internal_time) 
  {
    // Recurse toward the tree root on the infection tree, 
    // adding lineage breadcrumbs until we hit a preexisting genetic lineage 
    // or we hit the root.
    while(at(index).get_lineage() == NO_LINEAGE){
      at(index).set_lineage(new_leaf);
      if(at(index).mother() == MOTHER_EVE) break; 
      index = at(index).mother();
    }

    // Note the time of the internal node (or root)
    internal_time = at(index).time();

    // In the case of breaking a branch on the gene tree,
    // update breadcrumbs on infection tree to refer to new the internal node
    if(at(index).mother() != MOTHER_EVE) {
      gnode_index write_over_me = at(index).get_lineage();
      gnode_index internal_index = new_leaf + 1;
      while(at(index).get_lineage() == write_over_me){
	at(index).set_lineage(internal_index);
	if(at(index).mother() == MOTHER_EVE) break; 
	index = at(index).mother();
      }
      return write_over_me;
    }
    return -1;
  }

  /*
    INPUT: (1) Index of new leaf in gtree
           (2) Time of the new leaf the gene tree
	   (3) Index of the top node of the split branch (or -1 if no split)
           (4) Time of the new internal node (or existing root) in the gene tree
    MODIFIES: Adds a new branch to the gene tree. There are three cases:
                (1) Split an existing branch
		(2) First branch to be added to the tree
		(3) Branch to attach directly to the root
	      Note that (2) and (3) happen only once each.
	      This function updates everything on the gtree to properly add the 
	      new branch except the evolutionary time and the lmatrices. 
    OUTPUT: none
  */
  void add_gtree_branch(gnode_index new_leaf, 
			double leaf_time, 
			gnode_index split_branch_top,
			double internal_time)
  {
    if(split_branch_top != -1){
      
      //CASE 1: We need to break an existing branch on the gene tree
      //        in order to add the new branch. This means adding a 
      //        internal node as well as a new leaf.
      gnode_index internal_index = new_leaf + 1;
      gnode_index internal_mom = _gtree.at(split_branch_top).mother();
      // New leaf
      _gtree.push_back(gnode(internal_index, leaf_time)); 
      // New internal node
      _gtree.push_back(gnode(internal_mom, internal_time, new_leaf, 
			     split_branch_top)); 

      // Lastly, rewire connections in the broken branch
      if(_gtree.at(internal_mom).daughter_one() == split_branch_top) {
	_gtree.at(internal_mom).daughter_one(internal_index);			
      } else {
	_gtree.at(internal_mom).daughter_two(internal_index);			
      }
      _gtree.at(split_branch_top).mother(internal_index);

    } else if (split_branch_top == -1 && _gtree.size() == 1){

      // CASE 2: This is the first sample to be added to the tree
      // NOTE: First leaf from the root is always placed at daughter_one
      // First leaf
      _gtree.push_back(gnode(0, leaf_time)); 
      _gtree.at(0).daughter_one(new_leaf);

    } else if (split_branch_top == -1 && _gtree.size() > 0) {

      //CASE 3: There is the second sample to attach to the gene 
      // NOTE: Second leaf from the root is always placed at daughter_two,
      //       and the true root is always at index 1 in the gtree
      _gtree.push_back(gnode(0, leaf_time)); 
      _gtree.at(0).daughter_two(new_leaf);
    }
  }

  /*
    INPUT: (1) A gene tree
           (2) The parameters
	   (3) The index of the new leaf in the gene tree
    MODIFIES: Adds a random evolutionary branch length for the branch subtending
              the new leaf. This branch length is drawn from a gamma distribution 
	      with a mean equal to the branch length in calendar time. 
    OUTPUT: The length of the random branch
  */
  double gamma_branch(gsl_rng * rngptr, 
		      vector<gnode> & gtree,
		      map<string, double> & params, 
		      gnode_index new_leaf)
  {
    // Add gamma noise to the calendar time branch length
    gnode_index mom = gtree.at(new_leaf).mother();
    double time_branch = gtree.at(new_leaf).time() - gtree.at(mom).time(); 

    if(params["relax_branch"] == 0){
      
      // In the case of a degenerate relaxed clock, simply write branch length 
      // in calendar time to the evolutionary branch length slot
      time_branch += params["fixed_stem"];
      gtree.at(new_leaf).evo_branch_length(time_branch);
      return time_branch;
      
    } else {
      
      // Otherwise scale the calender branch length by multiplying by a gamma
      // random variable with mean one
      double evo_branch;

      if(time_branch == 0) {
	evo_branch = 0;
      } else {
	evo_branch = gsl_ran_gamma(rngptr,
				   time_branch / params["relax_branch"],
				   params["relax_branch"]);
      }
      evo_branch += params["fixed_stem"];
      gtree.at(new_leaf).evo_branch_length(evo_branch);
      return evo_branch;
    }
  }

  /*
    INPUT: (1) A gene tree
           (2) The sequences
           (3) The index on the gene tree of node with the lmatrix to recompute
    MODIFIES: Recomputes the lmatrix at fix_me. See below for the two cases.
    OUTPUT: none
  */
  void recompute_lmatrix(vector<gnode> & gtree, 
			 const vector<string> & sequences, 
			 gnode_index fix_me) 
  {
    // Extract the branch length subtending fix_me
    double evo_branch = gtree.at(fix_me).evo_branch_length();
    
    Lmatrix ell;
    if(gtree.at(fix_me).daughter_two() == NO_DAUGHTER) {
      // CASE 1: fix_me is TERMINAL, therefore we need to reattach its sequence
      // in order to recompute its lmatrix
      int seq_index = gtree.at(fix_me).daughter_one();
      ell = _model.sequence(sequences.at(seq_index));
      _model.backward_action(ell, evo_branch);
    } else {
      // CASE2 : fix_me is INTERNAL, therefore we need to merge its two daughter 
      // lmatrices in order to recompute its lmatrix
      node_index daughter_one = gtree.at(fix_me).daughter_one();
      node_index daughter_two = gtree.at(fix_me).daughter_two();
      ell = *(gtree.at(daughter_one).ellmatrix());
      ell *= *(gtree.at(daughter_two).ellmatrix());
      _model.backward_action(ell, evo_branch);
    }
    // Write the recomputed lmatrix to fix_me
    gtree.at(fix_me).ellmatrix(ell);
  }

  /*
    INPUT: (1) A gene tree
           (2) The parameters
	   (4) The index in the gene tree at the top of the split edge
    MODIFIES: Note that this function assumes calendar times are all correct and
              evolutionary branch lengths still need to be adjusted. 
	      This function divides the evolutionary time of the split edge to
	      the two new edges according to a beta distribution. 
    OUTPUT: none
  */
  void beta_bridge(gsl_rng * rngptr, 
		   vector<gnode> & gtree,
		   map<string, double> & params, 
		   gnode_index top) 
  {
    
    // Extract branch lengths in time on either side of the split as well as
    // the total evolutionary branch length to be divided
    gnode_index middle = gtree.at(top).mother();
    gnode_index root = gtree.at(middle).mother();
    double branch_one = gtree.at(top).time() - gtree.at(middle).time(); 
    double branch_two = gtree.at(middle).time() - gtree.at(root).time();
    double evo_branch = gtree.at(top).evo_branch_length(); 

    // Beta bridge to allocate branch lengths to the two new edges
    if(params["relax_branch"] == 0) {

      // In the case of a degenerate relaxed clock, simple write calendar time
      // branch lengths to the evolutionary branch length slots
      gtree.at(top).evo_branch_length(branch_one);
      gtree.at(middle).evo_branch_length(branch_two);

    } else {

      double evo_branch_one;
      double evo_branch_two;

      if(branch_one == 0 & branch_two == 0){

	evo_branch_one = 0;
	evo_branch_two = 0;

      } else if(branch_one == 0) {

	evo_branch_one = 0;
	evo_branch_two = evo_branch;

      } else if (branch_two == 0) {

	evo_branch_one = evo_branch;
	evo_branch_two = 0;

      } else {
	
	// Otherwise we have a nondegenerate beta bridge
	double proportion;
	proportion = gsl_ran_beta(rngptr,
				  branch_one / params["relax_branch"],
				  branch_two / params["relax_branch"]);    
	evo_branch_one = proportion * evo_branch;
	evo_branch_two = (1 - proportion) * evo_branch;
	
      }

      // Write evolutionary branch lengths
      gtree.at(top).evo_branch_length(evo_branch_one);
      gtree.at(middle).evo_branch_length(evo_branch_two);
    }
  }

  /*
    INPUT: (1) A vector of gtrees
           (2) A pointer to the start of an array of likelihoods with the same
	       length as the vector of gtrees
	   (3) The number of gtrees
    MODIFIES: Samples a single gtree from the vector of gtrees weighted by their
              likelihoods. Writes this sampled gtree to _gtree.
    OUTPUT: Returns the log of the average likelihood of the genetic trees held
            in gtrees.
  */
  double sample_gtree(gsl_rng * rngptr, 
		      const vector< vector<gnode> > & gtrees, 
		      double * logliks,
		      const int & num_samples)
  {
    // Scale log likelihoods by maximum, exponentiate
    double max_loglik = logliks[0];
    for (int i = 1; i < num_samples; i++)
      max_loglik = (logliks[i] > max_loglik) ? logliks[i] : max_loglik;
    for (int i = 0; i < num_samples; i++)
      logliks[i] = exp(logliks[i] - max_loglik);
    
    // Sample one genetic tree 
    for (int i = 1; i < num_samples; i++)
      logliks[i] += logliks[i-1];
    double du = logliks[num_samples - 1] / double(num_samples);
    double u = gsl_runif(rngptr, 0, logliks[num_samples - 1]);
    int sampled = 0;
    while (u > logliks[sampled]) sampled++;
    
    // Copy sampled tree into member variable _gtree
    _gtree = gtrees.at(sampled);
    
    // Update the cumulative likelihood
    _loglik += log(du) + max_loglik;

    // Return the average conditional log likelihood
    return log(du) + max_loglik;
  }

  /*
    INPUT: (1) The parameters
           (2) The sequences
	   (3) The index of the sequence to attach
	   (4) The index of the new leaf on the gtree
	   (5) The index of the top node of the split branch 
	       (or -1 if no split)
    MODIFIES: Either updates _gtree directly and then peels (in the case 
              of relax_branch == 0). Or, makes copies of the gene tree, 
	      tries num_branch_samples gamma distributed branch lengths, 
	      and performs a weighted sample to obtain a single gene tree 
	      to write over _gtree (in the case of relax_branch > 0).
    OUTPUT: Returns the conditional log likelihood of the particle
  */
  double relax_branch_lengths(gsl_rng * rngptr,
			      map<string, double> & params,
			      const vector<string> & sequences, 
			      int seq_index,
			      gnode_index new_leaf,
			      gnode_index split_branch_top) 
  {
    //Make copies of _gtree
    int num_samples = params["num_branch_samples"];
    if(params["relax_branch"] == 0 || seq_index == 0) num_samples = 1;
    vector< vector<gnode> > gtrees;
    gtrees.reserve(num_samples);
    for(int i = 0; i < num_samples; i++) gtrees.push_back(_gtree);
    double logliks[num_samples];      
    Lmatrix ell;
    double evo_branch;
    // CASE 1: We have broken a branch in the genetic tree and need to do a 
    //         beta bridge as well as sample branch lengths
    // CASE 2: We are connecting a new leaf directly to the root, and only 
    // 	 need to sample branch lengths
    for(int i = 0; i < num_samples; i++) {
      // CASE 1
      if(split_branch_top != -1) {
	beta_bridge(rngptr, gtrees.at(i), params, split_branch_top); 
	// Update the lmatrix at the top node of the split edge
	recompute_lmatrix(gtrees.at(i), sequences, split_branch_top);
      }

      // CASES 1 & 2
      evo_branch = gamma_branch(rngptr, gtrees.at(i), params, new_leaf); 

      // Attach the sequence at the new leaf
      ell = _model.sequence(sequences.at(seq_index));
      _model.backward_action(ell, evo_branch);
      gtrees.at(i).at(new_leaf).ellmatrix(ell);
      // Peel and store the likelihood
      logliks[i] = peel(gtrees.at(i), new_leaf);
    }
      
    // Sample a gene tree weighted by the likelihood
    if(num_samples == 1) {
      _gtree = gtrees.at(0);
      _loglik += logliks[0];
      return logliks[0];
    } else {
      return sample_gtree(rngptr, gtrees, logliks, num_samples);
    }
  }

  /*
    INPUT: (1) The index of the meristem to attach a sequence to
	   (2) The time of the sequence
	   (3) The parameters
	   (4) The index of the new leaf in the gene tree
    MODIFIES: The infection tree and the gene tree so as to allow for 
              either simulating a sequence or attaching a sequence.
	      Updates everything except evolutionary branch lengths 
	      and lmatrices. 
    OUTPUT: Either index of a new internal node in the gene tree or -1 (in
            the case of not splitting a branch)
  */
  gnode_index add_sequence_node(const meristem_index & sequenced_index, 
				double seq_time, 
				map<string,double> & params,
				const gnode_index new_leaf)
  {
    // Calculate the time to fix the leaf on both the infection tree and the 
    // genetic tree 
    double time_infected = at(sequenced_index).time();
    double new_leaf_time = 
      seq_time + params["stem_proportion"] * (seq_time - time_infected); 
    
    // Add leaf to represent the sequence on the infection tree
    node_index index = add_leaf(sequenced_index, seq_time, new_leaf_time);

    // Push sequence node to tree table if requested
    if(_make_tree_table) {
      node_index mom = at(index).mother();
      add_tree_table_row(index, mom, index, new_leaf_time);
    }

    // Add breadcrumbs for the new genetic lineage and write over other 
    // breadcrumbs if we have split a branch on the genetic tree
    double internal_time; // to be assigned inside update_breadcrumbs
    gnode_index split_branch_top = update_breadcrumbs(index, new_leaf, internal_time);
    
    // Add branch to the genetic tree, adjusting everything except evolutionary
    // branch length
    add_gtree_branch(new_leaf, 
		     new_leaf_time, 
		     split_branch_top, 
		     internal_time);

    // Return index of genetic node at the top of the split branch 
    // (or -1 in the case of no split branch)
    return split_branch_top;
  }
  
  /*
    INPUT: (1) The index of the meristem to attach a sequence to
           (2) The sequences
	   (3) The time of the sequence
	   (4) The index of the sequence
	   (5) The parameters
    MODIFIES: The infection tree and the gene tree. Details of each step 
              are documented for the series of functions that update each
	      structure. 
    OUTPUT: The conditional log likelihood of the attached sequence
  */
  double attach_sequence(gsl_rng * rngptr,
			 const meristem_index & sequenced_index, 
			 const vector<string> & sequences, 
			 const double seq_time, 
			 const int seq_index, 
			 map<string,double> & params)   
  {
    // Get the index of the new node to add to the gene tree
    gnode_index new_leaf = _gtree.size(); 

    // Add the new sequence node
    gnode_index split_branch_top = add_sequence_node(sequenced_index, 
						     seq_time, 
						     params,
						     new_leaf);
    
    // Add the sequence index to the new gene tree leaf
    _gtree.at(new_leaf).daughter_one(seq_index);

    // Try a number of branch lengths, sample a tree and return the log of the
    // average conditional likelihood
    return relax_branch_lengths(rngptr, params, sequences, seq_index, 
				new_leaf, split_branch_top);
  }

  /*
    INPUT: (1) A genetic tree
           (2) A node to start peeling from
    MODIFIES: Recurses down the passed gtree, updating lmatrices as it
              goes. Note that this function assumes that the lmatrix at 
	      the leaf node has already been passed through backward_action and
	      it assumes that internal lmatrices that do not lie directly on the 
	      path to the root have been properly updated.
    OUTPUT: The conditional log likelihood of the newly attached sequence
  */
  double peel (vector<gnode> & gtree, node_index me)
  {
    node_index mom = gtree.at(me).mother();
    node_index sis;
    if (gtree.at(mom).daughter_two() == me) {    
      sis = gtree.at(mom).daughter_one();
    } else {
      sis = gtree.at(mom).daughter_two();
    }
    Lmatrix *my_ell = gtree.at(me).ellmatrix();
    Lmatrix *moms_ell = gtree.at(mom).ellmatrix();   
    double new_loglik, old_loglik = 0; 
    if (is_groot(gtree, mom)) {
      old_loglik = moms_ell->loglik();
    }
    *moms_ell = *my_ell;
    if (sis != NO_LINEAGE){
      *moms_ell *= *(gtree.at(sis).ellmatrix());
    }  
    if (is_groot(gtree, mom)) { 
      _model.calc_prob(moms_ell);
      new_loglik = moms_ell->loglik();
      return new_loglik - old_loglik;
    } else {
      gnode_index grandma = gtree.at(mom).mother();
      double evo_branch_length = gtree.at(mom).evo_branch_length();
      _model.backward_action(*moms_ell, evo_branch_length);
      return peel(gtree,mom);
    }
  }
  
  /*
  INPUT: (1) An index of a meristem in the transmission tree
         (2) The number of loci to simulate
	 (3) The time of the sequence
	 (4) The parameters
  MODIFIES: Simulates a sequence conditional on all other sequences
            that have been simulated up to this point in time. 
	    Updates the transmission and genetic trees accordingly.
  OUTPUT: The simulated sequence
  */
  string simulate_sequence(gsl_rng * rngptr, 
			   meristem_index sequenced_index, 
			   int nlocus, 
			   double seq_time, 
			   map<string,double> & params) 
  {
    // Get the index of the new node to add to the gene tree
    gnode_index new_leaf = _gtree.size(); 

    // Add the new sequence node 
    // NOTE: This function adds nodes such that the calendar times
    //       are correct in the infection tree and the gene tree.
    //       (evolutionary branch lengths yet to be set)
    gnode_index split_branch_top = add_sequence_node(sequenced_index, 
						     seq_time, 
						     params,
						     new_leaf);
    // As in peeling, there are two cases:
    // CASE 1 (we have split a branch)
    if(split_branch_top != -1) {
      beta_bridge(rngptr, _gtree, params, split_branch_top); 
      bridge_lineages(rngptr, split_branch_top);
    }
    // CASES 1 & 2 (split branch or no split branch)
    double not_used = gamma_branch(rngptr, _gtree, params, new_leaf); 
    // Simulate the sequence at the new leaf  
    return _model.sequence(dress(rngptr, new_leaf, nlocus));
  }

//80////////////////////////////////////////////////////////////////////////////

/*
  INPUT: (1) The index of the leaf of the split branch
  MODIFIES: Simulates a sequence and adds it to the new internal node
  OUTPUT: none
 */
  void bridge_lineages (gsl_rng * rngptr, gnode_index leaf) {
  
  // Extract indices of new internal node and its mother
  gnode_index internal = _gtree.at(leaf).mother();
  gnode_index root = _gtree.at(internal).mother();

  // Get branchlengths to and from internal node
  double past_evo_branch = _gtree.at(internal).evo_branch_length();
  double future_evo_branch = _gtree.at(leaf).evo_branch_length();

  // Compute distribution for each base at the internal node given the states 
  // above and below
  Lmatrix *root_seq = _gtree.at(root).ellmatrix();	  
  Lmatrix *leaf_seq = _gtree.at(leaf).ellmatrix();	  
  Lmatrix *past(root_seq);
  Lmatrix *future(leaf_seq);
  _model.forward_action(past, past_evo_branch);
  _model.backward_action(*future, future_evo_branch);
  
  // Elementwise of multiplication past and future probabilities 
  // (write result over 'past' lmatrix)
  int nloci = past->nlocus;
  int nalleles = past->nallele;
  storage_type *p = past->data;
  storage_type *f = future->data;
  double *sum_array = new double[nloci];
  double *s = sum_array;
  int i, j;
  for (j = 0; j < nloci; j++) {
    double sum = 0;
    for (i = 0; i < nalleles; i++, p++, f++) {
      *p *= *f;	
      sum += *p;
    }
    *s++ = sum;
  }    
  
  // Normalize probabilites 
  for (j = 0, s = sum_array, p = past->data; j < nloci; j++) {
    for (i = 0; i < nalleles; i++, p++){
      *p /= *s;		
    }
    s++;  
  }  
  
  delete[] sum_array;
  
  // Simulate and place sequence
  _model.sim(rngptr, past);
  _gtree.at(internal).ellmatrix(*past);
}

/*
  INPUT: (1) An index of a node in the genetic tree
         (2) The number of loci to simulate
  MODIFIES: If a sequence does not already exist at the
            genetic node, this function simulates a sequence
	    and places it at the genetic node.
  OUTPUT: The simulated lmatrix
 */
  Lmatrix * dress(gsl_rng * rngptr, gnode_index me, int nlocus) {
    if (_gtree.at(me).has_data()) {
      return _gtree.at(me).ellmatrix();
    } else if (is_groot(me)) {
      Lmatrix *ell = _model.stat_ellmatrix(nlocus);
      _model.sim(rngptr, ell);
      _gtree.at(me).ellmatrix(*ell);
      return ell;
    } else {
      node_index mom = _gtree.at(me).mother();
      Lmatrix *ell = new Lmatrix(*(dress(rngptr, mom, nlocus)));
      double evo_branch = _gtree.at(me).evo_branch_length();
      _model.forward_action(ell, evo_branch);
      _model.sim(rngptr, ell);
      _gtree.at(me).ellmatrix(*ell);
      return ell;
    }
  }

  /*
  INPUT: (1) An index of a node in the genetic tree
         (2) The number of loci to simulate
  MODIFIES: If a sequence does not already exist at the
            genetic node, this function simulates a sequence
	    and places it at the genetic node.
  OUTPUT: The simulated lmatrix
 */
  void set_root_sequence(string root_sequence) {
    Lmatrix ell = _model.sequence(root_sequence);
    _gtree.at(0).ellmatrix(ell);
  }

/*
string tree::newick (void) {
  string s = newick(_roots.front());
  s += ";";
  return s;
}

string newick (node_index root) {
  stringstream ss;
  string s;
  double t;
  node_index mom, daughter[2];
  int d = 0;
  for (int i = 0; i < int(size()); i++) {
    if (at(i).mother() == root) 
      daughter[d++] = i;
  }
  mom = at(root).mother();
  t = (isnan(at(root).get_time())) ? _depth : at(root).get_time();
  t = (this->is_root(root)) ? t : t - at(mom).get_time();
  //if(t == 0) t = std::numeric_limits< double >::min();
  if (d == 0) {
    ss << root << ":" << t;
  } else {
    //cout << "daughter 1 of " << root << ":" << daughter[0] << ", time = " << t << endl;
    //cout << "daughter 2 of " << root << ":" << daughter[1] << ", time = " << t << endl;
    string s1, s2;
    s1 = newick(daughter[0]);
    s2 = newick(daughter[1]);
    s1 += ",";
    s1 += s2;
    ss << "(" << s1 << ")" << root << ":" << t;
  }
  ss >> s;
  return s;
}

string gnewick (node_index root = 1) {
    stringstream ss;
    string s;
    double t;
    node_index mom, daughter[2];
    int d = 0;
    for (int i = 0; i < int(_gtree.size()); i++) {
      if (_gtree.at(i).mother() == root) 
        daughter[d++] = i;
    }
    mom = _gtree.at(root).mother();
    t = _gtree.at(root).get_time();
    t = (_gtree.at(root).mother() == MOTHER_EVE) ? t : t - _gtree.at(mom).get_time();
    if (d==0) {
      ss << root << ":" << t;
    } else {
      if (d==1) {
    	return gnewick(daughter[0]);
      } else {
        string s1, s2;
        s1 = gnewick(daughter[0]);
        s2 = gnewick(daughter[1]);
        s1 += ",";
        s1 += s2;
        ss << "(" << s1 << ")" << root << ":" << t;
      }
    }
    ss >> s;
    return s;
  }

void tree::gnode_times(gnode_index me, doubvec& nodeTimes) {

  // In the case of a terminal node or a root with only one descendant
  // we do not save the node time as it is not internal and therefore   
  // not a coalescent event
  
  if(_gtree.at(me).daughter_two() != NO_DAUGHTER) {
    //Save the time of the current node
    nodeTimes.push_back(_gtree.at(me).get_time());
    //Pass the function down the tree
    gnode_times(_gtree.at(me).daughter_one(),nodeTimes);
    gnode_times(_gtree.at(me).daughter_two(),nodeTimes);
  }
  //In the rare case of a root with only one descendant, ensure the function does not stall
  if(_gtree.at(me).daughter_two() == NO_DAUGHTER && _gtree.at(me).mother() == MOTHER_EVE) {
    gnode_times(_gtree.at(me).daughter_one(),nodeTimes);  
  } 
}
  
vector<doubvec> tree::gnode_times() {

  vector<doubvec> result;
  int numTrees = _groots.size();
  
  if(_groots.size() == 0) cout << "No Roots!" << endl;
  
  result.reserve(numTrees);
  
  for (vector<node_index>::iterator j = _groots.begin(), end = _groots.end(); j != end; j++) {
    doubvec gnodeTimes;
    gnode_times(*j,gnodeTimes);
    result.push_back(gnodeTimes);
  }
  
  return result;
}

vector<double> tree::gtree_mrca_times() {
  vector<double> result;
  int numGnodes = _gtree.size();
  result.reserve(numGnodes);
  double mrcaTime;
  int mom;
  for(int i = 0; i != numGnodes; i++) {
    if(_gtree.at(i).daughter_one() == NO_DAUGHTER && _gtree.at(i).daughter_two() == NO_DAUGHTER) {
      mom = _gtree.at(i).mother();
      mrcaTime = _gtree.at(i).get_time() - _gtree.at(mom).get_time();
      result.push_back(mrcaTime);
    }
  }
  return result;
}
  */

void save_gtree(const char * outfile){
    ofstream out;
    out.open(outfile, std::ios_base::app);
    //print header
    out << "node ancestor time evobl\n";
    for(int i = 0; i < _gtree.size(); i ++){
      out << i << ' ' << _gtree.at(i).mother() << ' ' 
	  << _gtree.at(i).time() << ' ' 
	  << _gtree.at(i).evo_branch_length() << endl;
    }
    out.close();
}

void save_gtree(ofstream & outfile){
    //print header
    outfile << "node ancestor time evobl\n";
    for(int i = 0; i < _gtree.size(); i ++){
      outfile << i << ' ' << _gtree.at(i).mother() << ' ' 
	      << _gtree.at(i).time() << ' ' 
	      << _gtree.at(i).evo_branch_length() << endl;
    }
}

void save_transmission_tree(ofstream & outfile){
    //print header
    outfile << "node ancestor time glineage\n";
    for(int i = 0; i < size(); i ++){
      outfile << i << ' ' << at(i).mother() << ' ' 
	      << at(i).time() << ' ' 
	      << at(i).get_lineage() << endl;
    }
}

void save_transmission_tree(const char * outfile){
    ofstream out;
    out.open(outfile, std::ios_base::app);
    //print header
    out << "node ancestor time glineage\n";
    //print node info
    for(int i = 0; i < size(); i ++){
      out << i << ' ' << at(i).mother() << ' ' 
	  << at(i).time() << ' ' 
	  << at(i).get_lineage() << endl;
    }
    out.close();
}

};

#endif
