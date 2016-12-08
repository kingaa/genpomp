#include "gsl_rng.h"

/*
Returns a random integer from the generator r, with all values between min and max 
equally likely.
*/
unsigned long int
gsl_rng_get (gsl_rng * r)
{
  return gsl_rng_uniform_int(r, RNGSTREAM_MAX);
}

/*
Returns a random number in the interval (0,1).
*/
double
gsl_rng_uniform (gsl_rng * r)
{
  return RngStream_RandU01(*r);
}

/*
Returns a random number in the interval (0,1).
Note: when using L'ecuyer's RngStream rng this function should be the 
same as gsl_rng_uniform.
*/
double
gsl_rng_uniform_pos (gsl_rng * r)
{
  double x ;
  do
    {
      x = RngStream_RandU01(*r);
    }
  while (x == 0) ;

  return x ;
}

/*
Returns a random integer in the interval [0, n-1].
*/
unsigned long int
gsl_rng_uniform_int (gsl_rng * r, unsigned long int n)
{
  int max = RNGSTREAM_MAX;
  if (n > max || n == 0) 
    {
      perror("invalid n, either 0 or exceeds maximum value of generator");
    }
  return RngStream_RandInt (*r, 0, n - 1);
}
