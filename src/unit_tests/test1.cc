// -*- mode: C++; -*-

/*
TEST1

Constructs Tamura Nei molecular model 
Checks that the matrix is properly constructed

*/

#define MATHLIB_STANDALONE 
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <time.h>
#include "lmatrix.h"
#include "userUnif.h"
#include "substmodel.h"
#include "tree.h"
#include "abstractUsermodel.h"
#include "io.h"
#include "reduceSeqs.h"
#include "type_defs.h"

using namespace std;

int main ()
{
  //parameters that throw errors
  map<string,double> params;
  params["beta"] = 0.1;
  params["alphaR"] = 0.05;
  params["alphaY"] = 0.002;
  params["piA"] = 0.4;
  params["piG"] = 0.2;
  params["piT"] = 0.2;
  params["piC"] = 0.2;
  params["relax"] = 0.1;

  for(map<string, double>::const_iterator it = params.begin(); it != params.end(); ++it)
    {
      cout << it->first << " " << it->second << "\n";
    }
   
  //Throw in some values here just for initializing
  int numPatterns = 100;
  double patternCountArray[numPatterns];
  for(int i = 0; i < numPatterns; i++) patternCountArray[i] = 1;

  // Set molecular model of sequence evolution
  double Q[16];
  string molModelName = "tamuraNei";
  makeQmatrix(molModelName,params,Q);
  substModel mod(4, Q, "AGCT", numPatterns, patternCountArray, params["relax"]);
  
  cout << endl;
  cout << Q[0] << ' ' << Q[4] << ' ' << Q[8] << ' ' << Q[12] << endl;
  cout << Q[1] << ' ' << Q[5] << ' ' << Q[9] << ' ' << Q[13] << endl;
  cout << Q[2] << ' ' << Q[6] << ' ' << Q[10] << ' ' << Q[14] << endl;
  cout << Q[3] << ' ' << Q[7] << ' ' << Q[11] << ' ' << Q[15] << endl;
  cout << endl;

}
