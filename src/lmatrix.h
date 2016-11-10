// -*- mode: C++; -*-
#ifndef _LMATRIX_H_
#define _LMATRIX_H_

#include <cstring>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "type_defs.h"

using namespace std;

extern "C" {
  void merge_lineages(int n, int nlocus,
		      storage_type *my_data, double *my_charac,
		      const storage_type *her_data, double her_charac);
}

class Lmatrix {
public:
  int nallele, nlocus;
  double characteristic;
  storage_type *data;

  Lmatrix (void) : nallele(0), nlocus(0), characteristic(0), data(NULL) { }

  Lmatrix (int na, int nl, const storage_type *dat) : nallele(na), nlocus(nl),
    characteristic(0), data(NULL) { 
    if (na*nl > 0) {
      data = new storage_type[na*nl];
      memcpy(data,dat,na*nl*sizeof(storage_type));
    }
  }

  Lmatrix (const Lmatrix& ell) : nallele(ell.nallele), nlocus(ell.nlocus), 
    characteristic(ell.characteristic), data(NULL) { 
    if (ell.size() > 0) {
      data = new storage_type[ell.size()];
      memcpy(data,ell.data,ell.size()*sizeof(storage_type));
    }
  }

  Lmatrix& operator= (const Lmatrix& ell) { 
    if (this->size() != ell.size()) {
      this->nallele = ell.nallele;
      this->nlocus = ell.nlocus;
      this->characteristic = ell.characteristic;
      if (this->data != NULL) delete[] this->data;
      if (ell.size() > 0) {
    	this->data = new storage_type[ell.size()];
	memcpy(this->data,ell.data,ell.size()*sizeof(storage_type));
      } else {
    	this->data = NULL;
      }
    } else {
      this->nallele = ell.nallele;
      this->nlocus = ell.nlocus;
      this->characteristic = ell.characteristic;
      if (this->data == NULL) this->data = new storage_type[ell.size()];
      if (ell.size() > 0)
	memcpy(this->data,ell.data,ell.size()*sizeof(storage_type));
    }
    return *this;
  }

  ~Lmatrix (void) {
    if (data != NULL) delete[] data;
  }

  int size (void) const {
    return nallele*nlocus;
  }

  storage_type loglik (void) const {
    return characteristic;
  }

  Lmatrix& operator*= (const Lmatrix& ell) {
    
    if (this->nallele != ell.nallele) {
      cout << "nalleles not the same:" << endl;
      cout << this->nallele << endl;
      cout << ell.nallele << endl;
    }
    
    if (this->nlocus != ell.nlocus) {
      cout << "nlocus not the same:" << endl;
      cout << this->nlocus << endl;
      cout << ell.nlocus << endl;
    }
    
    if (this->nallele != ell.nallele || this->nlocus != ell.nlocus)
      throw(runtime_error("incommensurate Lmatrix dimensions"));
    merge_lineages(this->nallele, this->nlocus,
		   this->data, &(this->characteristic),
		   ell.data, ell.characteristic);
    return *this;
  }

  friend ostream& operator<< (ostream& o, Lmatrix& ell) {
    o << "characteristic: " << ell.characteristic << endl;
    for (int i = 0; i < ell.nallele; i++) {
      o << ell.data[i];
      for (int j = 1; j < ell.nlocus; j++)
        o << "," << ell.data[i+ell.nallele*j];
      o << endl;
    }
    return o;
  }

};

#endif
