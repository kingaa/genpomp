#include "RngStream.h"
#include "gsl_rng.h"

extern "C" {
  void destroyStreams (void);
  gsl_rng * createStreams (int n, int seed_in);
  void save_rng_state(int seed_index);
  void set_rng_state(int seed_index);
  void destroy_states(void);
  double gsl_runif(gsl_rng *stream, double a, double b);
}

