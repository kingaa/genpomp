extern "C" {
  void destroyStreams (void);
  void createStreams (int n, int nthreads, Int32 seed_in);
  void save_rng_state(int seed_index);
  void set_rng_state(int seed_index);
  void create_initial_rng_states(int n_rng_states, Int32 seed_in);
  void destroy_states(void);
}

