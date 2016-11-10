#include "reduceSeqs.h"

// Function to compare strings for sorting purposes
bool stringCompare( const string &left, const string &right ){
  for( string::const_iterator lit = left.begin(), rit = right.begin(); lit != left.end() && rit != right.end(); ++lit, ++rit )
    if( tolower( *lit ) < tolower( *rit ) )
      return true;
    else if( tolower( *lit ) > tolower( *rit ) )
      return false;
  if( left.size() < right.size() )
    return true;
  return false;
}

// Function to transpose sequences
// Takes a vector of sequences (note: assumes that all sequences are the same length)
// Returns a vector of transposed sequences
// Assumes all sequences are the same length 
vector<string> transpose_seqs( const vector<string> &seqs) {
  vector<string> transposed;
  string t;
  for(int i = 0; i < seqs.at(0).length(); i++) {
    t.clear();
    for(int j = 0; j < seqs.size(); j++) {
      t += seqs.at(j).at(i);
    }
    transposed.push_back(t);
  }
  return(transposed);
}

// Function to remove the pattern of all missing (if it is present)
void remove_invariant_pattern(vector<string> & seq_patterns, vector<double> & counts, char invariant){
  int seqLength = seq_patterns.at(0).length();
  string useless (seqLength, invariant);
  for(int i = 0; i < seq_patterns.size(); i++){
    if(seq_patterns.at(i) == useless) {
      seq_patterns.erase(seq_patterns.begin() + i);
      counts.erase(counts.begin() + i);
    }
  }
}

// Function that writes over a vector of sequences (including NAs) with their reduced
// representation and then returns a vector of counts of how many times each column
// pattern is represented in the full sequence

vector<double> reduce_sequences(vector<string> & fullSeqs, bool removeInvariant){
  //First extract the sequences from the combined vector of sequences and NAs
  vector<string> trueSeqs;
  vector<int> trueSeqIndices;
  for(int i = 0; i < fullSeqs.size(); i++){
    if(fullSeqs.at(i).length() > 2){ // because NA has two characters
      trueSeqs.push_back(fullSeqs.at(i));
      trueSeqIndices.push_back(i);
    }
  }
  
  // Handle the case of no sequences
  if(trueSeqs.size() == 0){
    std::cout << "No sequences!" << std::endl;
    vector<double> one_count;
    one_count.push_back(1);
    return(one_count);
  }

  //Next transpose sequences and sort to allow for tallying column patterns
  vector<string> transposed = transpose_seqs(trueSeqs);
  sort( transposed.begin(), transposed.end(), stringCompare );

  //Loop through the sorted transposed sequences, tallying the count of each 
  vector<string> unique;
  vector<double> counts;
  double counter = 1;
  unique.push_back(transposed.at(0));
  for(int j = 1; j < transposed.size(); j++) {
    if(unique.back().compare(transposed.at(j)) == 0) {
	counter += 1;
      } else {
	unique.push_back(transposed.at(j));
	counts.push_back(counter);
	counter = 1;
    }
    if(j == transposed.size() -1) counts.push_back(counter);
  }    

  //Remove missing pattern if it exists (as this pattern contains no information)
  remove_invariant_pattern(unique, counts,'N');
  remove_invariant_pattern(unique, counts,'-');

  //If asked to, remove invariant sites
  if(removeInvariant){
    remove_invariant_pattern(unique, counts,'A');
    remove_invariant_pattern(unique, counts,'G');
    remove_invariant_pattern(unique, counts,'C');
    remove_invariant_pattern(unique, counts,'T');
  }

  //Transpose back, write over full sequences with reduced
  vector<string> reduced = transpose_seqs(unique);
  int j = 0;
  for(int i = 0; i < fullSeqs.size(); i++){
    if(i == trueSeqIndices.at(j)){
      cout << i << ' ' << reduced.at(j) << endl;
      fullSeqs.at(i) = reduced.at(j);
      j++;
    }
    if(j == reduced.size()) break;
  }
  cout << endl;

  for(int i = 0; i < counts.size(); i++) cout << counts.at(i) << ' ';
  cout << endl;

  //Finally, return pattern counts
  return(counts);
}
