#ifndef PADDED_MAP_H
#define PADDED_MAP_H

#ifdef  __cplusplus
extern "C" {
#endif

  typedef struct padded_map {
    int stream_index;
    double pad[16];
  } padded_map;
  
  
#ifdef  __cplusplus
}
#endif

#endif  /* NODE_H */  
