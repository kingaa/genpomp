#include "substmodel.h"
#include "type_defs.h"

void jukesCantor (double alpha, double * Q) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Q[(i*4) + j] = (i == j) ? -3*alpha : alpha;
    }
  }
}

void kimura (double alpha, double beta, double * Q) {
  // reduces to jukesCantor for
  // alpha = beta
  Q[0] = Q[5] = Q[10] = Q[15] = -alpha-2*beta;  //diagonals  
  Q[1] = Q[4] = Q[11] = Q[14] = alpha;
  Q[2] = Q[3] = Q[6] = Q[7] = beta;
  Q[8] = Q[9] = Q[12] = Q[13] = beta;
}

void F81 (double piA, double piG, double piC, double piT, double * Q) {
  Q[0] = -(piC + piA + piG);
  Q[1] = piC;  
  Q[2] = piA;  
  Q[3] = piG;
  
  Q[4] = piT;  
  Q[5] = -(piT + piA + piG);
  Q[6] = piA;  
  Q[7] = piG;  

  Q[8] = piT;  
  Q[9] = piC;  
  Q[10] =-(piT + piC + piG);
  Q[11] = piG;  

  Q[12] = piT;  
  Q[13] = piC;  
  Q[14] = piG;  
  Q[15] = -(piT + piC + piG);
}

void tamuraNei (double alphaR, double alphaY, double beta, double piA, double piG, double piC, double piT, double * Q) {
  //Rates of leaving T 
  Q[1] = alphaY*piC;  
  Q[2] = beta*piA;  
  Q[3] = beta*piG;
  Q[0] = -Q[1]- Q[2]- Q[3];

  //Rates of leaving C 
  Q[4] = alphaY*piT;  
  Q[6] = beta*piA;  
  Q[7] = beta*piG;  
  Q[5] = -Q[4] -Q[6] -Q[7];

  //Rates of leaving A 
  Q[8] = beta*piT;  
  Q[9] = beta*piC;  
  Q[11] = alphaR*piG;  
  Q[10] = -Q[8]- Q[9] -Q[11];

  //Rates of leaving G 
  Q[12] = beta*piT;  
  Q[13] = beta*piC;  
  Q[14] = alphaR*piA;  
  Q[15] = -Q[12]- Q[13]- Q[14];
}

void makeQmatrix(string molModel, map<string,double> &params, double * Q){
  int err = -1;
  if(molModel == "jukesCantor"){ jukesCantor(params["alpha"],Q); err++;}
  if(molModel == "kimura"){ kimura(params["alpha"],params["beta"],Q); err++;}
  if(molModel == "F81"){ F81(params["piA"],params["piG"],
			     params["piC"],params["piT"],Q); err++;}
  if(molModel == "tamuraNei"){ tamuraNei(params["alphaR"], params["alphaY"],
			       params["beta"],params["piA"],params["piG"],
			       params["piC"],params["piT"],Q); err++;}
  if(err != 0) cout << "Error: molecular model not properly named" << endl;
}

/*
NumericMatrix jukesCantor (double alpha) {
  NumericMatrix Q(4,4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Q(i,j) = (i == j) ? -3*alpha : alpha;
    }
  }
  return Q;
}

NumericMatrix kimura (double alpha, double beta) {
  NumericMatrix Q(4,4);
  // reduces to jukesCantor for
  // alpha = beta
  Q(0,0) = Q(1,1) = Q(2,2) = Q(3,3) = -alpha-2*beta;
  Q(0,1) = Q(1,0) = Q(3,2) = Q(2,3) = alpha;
  Q(0,2) = Q(0,3) = Q(1,2) = Q(1,3) = beta;
  Q(2,0) = Q(3,0) = Q(2,1) = Q(3,1) = beta;
  return Q;
}

NumericMatrix tamuraNei (double alphaR, double alphaY, double beta,
			 const double *pi) {
  // reduces to kimura for 
  // alphaR = alphaL = (alpha-beta)/2
  // pi[i] = 1/4 for all i
  NumericMatrix Q(4,4);
  double Pi[4];
  double sum = 0;
  for (int i = 0; i < 4; i++) sum += pi[i];
  for (int i = 0; i < 4; i++) Pi[i] = pi[i]/sum;
  double pr = Pi[0] + Pi[1];
  double py = Pi[2] + Pi[3];

  Q(0,1) = (alphaR/pr + beta)*Pi[0];
  Q(1,0) = (alphaR/pr + beta)*Pi[1];
  Q(0,2) = Q(0,3) = beta*Pi[0];
  Q(1,2) = Q(1,3) = beta*Pi[1];
  Q(2,0) = Q(2,1) = beta*Pi[2];
  Q(3,0) = Q(3,1) = beta*Pi[3];
  Q(2,3) = (alphaY/py + beta)*Pi[2];
  Q(3,2) = (alphaY/py + beta)*Pi[3];
  Q(0,0) = -Q(1,0)-Q(2,0)-Q(3,0);
  Q(1,1) = -Q(0,1)-Q(2,1)-Q(3,1);
  Q(2,2) = -Q(0,2)-Q(1,2)-Q(3,2);
  Q(3,3) = -Q(0,3)-Q(1,3)-Q(2,3);

  return Q;
}
*/
