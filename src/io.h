#ifndef _MY_IO_H_
#define __MY_IO_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "type_defs.h"

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
double StringToFloat_Type ( const string &Text);
void readParams( map<string,double> &paramMap, char *fileName);
void read_random_walk_params( map<string,double> &paramMap, vector<string> &transforms, vector<string> &perturbations, char *fileName);
void readSeqs( vector<double> &times, vector<string> &seqs, const char *seqFile);
void printStatus(ofstream &outfile, bool start, string functionName, int iter);

#endif

