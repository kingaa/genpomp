#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rdefines.h>
#include "RngStream.h"
#include "padded_map.h"


//double *u;        // array to hold one double for each RNG (note: unused in standalone version)

int nStreams = 0;
RngStream *rngstreams; // independent RNGs
//int *thread_stream_map = 0; 
padded_map *thread_stream_map = 0;


void destroyStreams (void) {
  int i;
  if (nStreams > 0) {
    for (i = 0;  i < nStreams;  i++) RngStream_DeleteStream(&rngstreams[i]);    
    free(rngstreams);
  }
  if (thread_stream_map != 0) free(thread_stream_map);
  nStreams = 0;
}

void createStreams (int n, int nthreads, Int32 seed_in) {
  int i;
  unsigned long lecseed[6];
  destroyStreams();
  if (n > 0) {
    nStreams = n;
    printf("Using %d RNG streams\n", n);
  } else if (nStreams == 0) {
    nStreams = 1;
    printf("Warning: using only one RNG\n");
  }

  for (i = 0; i < 6; i++) {
    lecseed[i] = (unsigned long) seed_in;
  }

  rngstreams = (RngStream *) malloc(nStreams*sizeof(RngStream *));
  //thread_stream_map = (int *) malloc(nthreads*sizeof(int));
  //for (i = 0; i < nthreads; i++) thread_stream_map[i] = 0;
  thread_stream_map = (padded_map *) malloc(nthreads*sizeof(padded_map));
  for (i = 0; i < nthreads; i++) thread_stream_map[i].stream_index = 0;

  int fail = RngStream_SetPackageSeed(lecseed);
  if (fail) printf("Warning: failure in RngStream_SetPackageSeed\n");
  for (i = 0;  i < nStreams;  i++) 
    rngstreams[i] = RngStream_CreateStream("");  
}

double unif_rand() 
{ 
  //return RngStream_RandU01(rngstreams[thread_stream_map[omp_get_thread_num()]]);
  return RngStream_RandU01(rngstreams[thread_stream_map[omp_get_thread_num()].stream_index]);
} 
