#ifndef _REDUCESEQS_H_
#define __REDUCESEQS_H_


#include <fstream>
#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "type_defs.h"

using namespace std;

bool stringCompare(const string &left, const string &right);
vector<string> transpose_seqs( const vector<string> &seqs);
void remove_invariant_pattern(vector<string> & seq_patterns, vector<double> & counts, char invariant);
vector<double> reduce_sequences(vector<string> & fullSeqs, bool removeInvariant);

#endif
