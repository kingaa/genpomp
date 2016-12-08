// -*- mode: C++; -*-
#ifndef _SUBST_MODEL_H_
#define _SUBST_MODEL_H_

#define MATHLIB_STANDALONE 


#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <cstring>
#include <string>
#include <stdexcept>
#include <map>
#include "lmatrix.h"
#include "type_defs.h"
#include "gsl_rng.h"
#include "gsl_randist.h"
#include "userUnif.h"

typedef std::map<char,int> allele_map;

using namespace std;

#define NO_DATA (-1)

extern "C" {
  void init_sm (int n, const double *qmatrix, 
		double *vals, 
		double *lv, double *rv, 
		double *bl, double *br, 
		double *sd,
		int *info, int *err);
  void molmodel_action (int n, int nlocus, double t, storage_type *ell,
			double sigma2,
			const double *vals, 
			const double *lv, const double *rv);
  void statdist_action (int nallele, int nlocus, storage_type *ell, 
			double *charac, const double *sd, const double *pattern_count);
}

class substModel {

private:
  int _nallele;
  string _astring;
  allele_map _amap;
  double *_qmatrix;  
  double *_statdist;
  double *_eigenvals;
  double *_leigenvec;
  double *_reigenvec;
  double *_reig_trans;
  double *_leig_trans;
  int _npatterns;         // the number of distinct locus patterns
  double *_pattern_count; // pointer to array with the number of times each locus pattern is represented in the set of full sequences
  double _relax;

  void _initialize (const double *qmatrix) {
    int info = 0, err = 0;

    init_sm(_nallele,_qmatrix,
	    _eigenvals,
	    _leigenvec,_reigenvec,
	    _reig_trans,_leig_trans,
	    _statdist,
	    &info,&err);
    if (info != 0)
      throw(runtime_error("error in eigendecomposition"));
    if (err == 1)
      throw(runtime_error("Q matrix is not stochastic"));
    if (err == 2)
      throw(runtime_error("Q matrix is not reversible"));
    if (err != 0)
      throw(runtime_error("bad Q matrix"));
  }

  void _allocate (int nallele, const double *qmatrix, 
		  const char *mstring, double relax) {
    _nallele = nallele;
    _qmatrix = new double[nallele*nallele];
    _statdist = new double[nallele];
    _eigenvals = new double[nallele];
    _leigenvec = new double[nallele*nallele];
    _reigenvec = new double[nallele*nallele];
    _reig_trans = new double[nallele*nallele];
    _leig_trans = new double[nallele*nallele];
    memcpy(_qmatrix,qmatrix,nallele*nallele*sizeof(double));
    _initialize(_qmatrix);
    for (int i = 0; i < nallele; i++) _amap[mstring[i]] = i;
    _astring = mstring;
    _relax = relax;
  }

  void _destroy (void) {
    delete[] _qmatrix;
    delete[] _statdist;
    delete[] _eigenvals;
    delete[] _leigenvec;
    delete[] _reigenvec;
    delete[] _reig_trans;
    delete[] _leig_trans;
  }

public:

  substModel (void) : _nallele(0), _qmatrix(NULL), _statdist(NULL),
    _eigenvals(NULL), _leigenvec(NULL), _reigenvec(NULL), 
		      _reig_trans(NULL), _leig_trans(NULL), _npatterns(0), _pattern_count(NULL), _relax(0) { 
    }
  
  substModel (int nallele, const double *qmatrix, const char *mstring, int npatterns, double *pattern_count, 
	      double relax = 0) {
    _npatterns = npatterns;
    _pattern_count = pattern_count;
    _allocate(nallele,qmatrix,mstring,relax);
  }

  substModel (const substModel& x) : _astring(x._astring), 
				     _amap(x._amap), _relax(x._relax), _qmatrix(NULL) {
    _nallele = x.nallele();
    _npatterns = x._npatterns;
    _pattern_count = x._pattern_count;

    if (_qmatrix != NULL) _destroy();
    _qmatrix = new double[_nallele*_nallele];
    _statdist = new double[_nallele];
    _eigenvals = new double[_nallele];
    _leigenvec = new double[_nallele*_nallele];
    _reigenvec = new double[_nallele*_nallele];
    _reig_trans = new double[_nallele*_nallele];
    _leig_trans = new double[_nallele*_nallele];
    memcpy(_qmatrix,x._qmatrix,_nallele*_nallele*sizeof(double));
    memcpy(_leigenvec,x._leigenvec,_nallele*_nallele*sizeof(double));
    memcpy(_reigenvec,x._reigenvec,_nallele*_nallele*sizeof(double));
    memcpy(_reig_trans,x._reig_trans,_nallele*_nallele*sizeof(double));
    memcpy(_leig_trans,x._leig_trans,_nallele*_nallele*sizeof(double));
    memcpy(_eigenvals,x._eigenvals,_nallele*sizeof(double));
    memcpy(_statdist,x._statdist,_nallele*sizeof(double));
  }

  substModel& operator= (const substModel& x) {
    _astring = x._astring;
     _amap = x._amap;
    _relax = x._relax;

    _npatterns = x._npatterns;
    _pattern_count = x._pattern_count;

    if (_nallele != x.nallele() || _qmatrix == NULL) {
      _nallele = x.nallele();
      if (_qmatrix != NULL) _destroy();
      _qmatrix = new double[_nallele*_nallele];
      _statdist = new double[_nallele];
      _eigenvals = new double[_nallele];
      _leigenvec = new double[_nallele*_nallele];
      _reigenvec = new double[_nallele*_nallele];
      _reig_trans = new double[_nallele*_nallele];
      _leig_trans = new double[_nallele*_nallele];
    }
    memcpy(_qmatrix,x._qmatrix,_nallele*_nallele*sizeof(double));
    memcpy(_leigenvec,x._leigenvec,_nallele*_nallele*sizeof(double));
    memcpy(_reigenvec,x._reigenvec,_nallele*_nallele*sizeof(double));
    memcpy(_reig_trans,x._reig_trans,_nallele*_nallele*sizeof(double));
    memcpy(_leig_trans,x._leig_trans,_nallele*_nallele*sizeof(double));
    memcpy(_eigenvals,x._eigenvals,_nallele*sizeof(double));
    memcpy(_statdist,x._statdist,_nallele*sizeof(double));
    return *this;
  }

  ~substModel (void) { 
    _destroy();
  }

  size_t size (void) const {
    int n = _nallele;
    return n*n;
  }

  size_t nallele (void) const {
    return _nallele;
  }

  int amap (const char c) {
    allele_map::iterator i = _amap.find(c);
    if (i==_amap.end()) return NO_DATA;
    else return i->second;
  }


  Lmatrix *stat_ellmatrix (int nlocus) {
    storage_type *sd = new storage_type[_nallele*nlocus];
    int i;
    for (i = 0; i < nlocus; i++)
      memcpy(&sd[i*_nallele],_statdist,_nallele*sizeof(storage_type));
    Lmatrix *ell = new Lmatrix(_nallele,nlocus,sd);
    delete[] sd;
    return ell;
  }

  /*
  void normalize (Lmatrix *ell) {
    int i, j;
    double sum, *ellp;
    for (j = 0, ellp = ell->data; j < ell->nlocus; j++) {
      for (i = 0, sum = 0; i < ell->nallele; i++)
	sum += ellp[i]*_statdist[i];
      for (i = 0; i < ell->nallele; i++, ellp++)
	*ellp *= _statdist[i]/sum;
    }
  }
  */

  void calc_prob (Lmatrix *ell) {
    if (_nallele != ell->nallele){
      cout << "Substmodel _nallele = " << _nallele << endl;
      cout << "Lmatrix _nallele = " << ell-> nallele << endl;
      throw(runtime_error("calc_prob: substModel and Lmatrix differ in number of alleles:"));
    }
    statdist_action(ell->nallele,ell->nlocus,ell->data,
		    &(ell->characteristic),_statdist,_pattern_count);
  }

  void forward_action (Lmatrix *ell, double t) {
    int nlocus = ell->nlocus;
    if (_nallele != ell->nallele)
      throw(runtime_error("forward_action: substModel and Lmatrix differ in number of alleles"));
    molmodel_action(_nallele, nlocus, t, ell->data, _relax,
		    _eigenvals, _leigenvec, _reigenvec);
  }

  void backward_action (Lmatrix& ell, double t) {
    int nlocus = ell.nlocus;
    if (_nallele != ell.nallele)
      throw(runtime_error("backward action: substModel and Lmatrix differ in number of alleles"));
    molmodel_action(_nallele, nlocus, t, ell.data, _relax,
		    _eigenvals, _reig_trans, _leig_trans);
  }

  Lmatrix sequence (const string& s) {
    int nlocus = s.length();
    storage_type *dat = new storage_type[_nallele*nlocus];
    storage_type *dp, *de;
    int i, j;
    de = dat + _nallele*nlocus;
    for (dp = dat; dp != de; dp++) *dp = 0;
    for (j = 0, dp = dat; j < nlocus; j++, dp += _nallele) {
      i = this->amap(s[j]);
      if (i != NO_DATA) 
        dp[i] = 1;
      else 
        for (i = 0; i < _nallele; i++) dp[i] = 1;
    }
    Lmatrix ell(_nallele,nlocus,dat);
    delete[] dat;
    return ell;
  }

  string sequence (Lmatrix *ell) {
    int n = ell->nallele;
    int nlocus = ell->nlocus;
    storage_type *ellp = ell->data;
    string seq(nlocus,'-');
    int i, j;
    for (j = 0; j < nlocus; j++) {
      for (i = 0; i < n; i++, ellp++) {
    if (*ellp == 1) seq[j] = _astring[i];
      }
    }
    return seq;
  }

///////////////////////////////////////////////
/// Was this correctly converted to cpp? //////
///////////////////////////////////////////////

  void sim (gsl_rng * rngptr, Lmatrix *ell) {
    int n = ell->nallele;
    int nlocus = ell->nlocus;
    storage_type c, *ellp = ell->data;
    double u;
    int i, j, k;
    for (j = 0; j < nlocus; j++) {
      u = gsl_runif(rngptr, 0, 1);
      c = ellp[0]; k = 0;
      while (c < u) c += ellp[++k];
      for (i = 0; i < n; i++)
	*(ellp++) = (i == k) ? 1 : 0;
    }
  }

};

void jukesCantor (double alpha, double * Q);
void kimura (double alpha, double beta, double * Q);
void F81 (double piA, double piG, double piC, double piT, double * Q);
void tamuraNei (double alphaR, double alphaY, double beta, double piA, double piG, double piC, double piT, double * Q);
void makeQmatrix(string molModel, map<string,double> &params, double * Q);

#endif
