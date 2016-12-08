#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "RngStream.h"
#include "gsl_rng.h"

/*
typedef struct paddedRngptr {
  RngStream rng; // 8 bytes
  double pad[15]; // 120 bytes 
} paddedRngptr;  // 128 bytes total = size of two cache lines
//double *u;        // array to hold one double for each RNG (note: unused in standalone version)
//paddedRngptr * rngstreams; // pointer to array of padded pointers to eliminate false sharing
*/

int nStreams = 0;
gsl_rng *rngstreams; // independent RNGs

void destroyStreams (void) {
  int i;
  if (nStreams > 0) {
    for (i = 0;  i < nStreams;  i++) RngStream_DeleteStream(&rngstreams[i]);    
    free(rngstreams);
  }
  nStreams = 0;
}

gsl_rng * createStreams (int n, int seed_in) {
  int i;
  unsigned long lecseed[6];
  destroyStreams();
  if (n > 0) {
    nStreams = n;
    printf("Using %d random number streams\n", n);
  } else if (nStreams == 0) {
    nStreams = 1;
    printf("Warning: using only one random number stream\n");
  }

  for (i = 0; i < 6; i++) {
    lecseed[i] = (unsigned long) seed_in;
  }

  rngstreams = (RngStream *) malloc(nStreams*sizeof(RngStream *));

  int fail = RngStream_SetPackageSeed(lecseed);
  if (fail) printf("Warning: failure in RngStream_SetPackageSeed\n");
  for (i = 0;  i < nStreams;  i++) 
    rngstreams[i] = RngStream_CreateStream("");  
  return rngstreams;
}

double gsl_runif(gsl_rng *stream, double a, double b)
{
  if (!isfinite(a) || !isfinite(b) || b < a)  NAN;

  if (a == b)
    return a;
  else {
    double u;
    do {u = RngStream_RandU01(*stream);} while (u <= 0 || u >= 1);
    return a + (b - a) * u;
  }
}

















