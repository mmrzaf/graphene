#ifndef GRAPHENE_UTILS_H
#define GRAPHENE_UTILS_H

// Random number generation (deterministic)
void utils_seed_rng(unsigned int seed);
int utils_rand_int(int min, int max);
double utils_rand_double(void);

// Array utilities
void utils_shuffle_int_array(int* arr, int n);
int* utils_load_seedset(const char* path, int* out_count);
void utils_save_seedset(const char* path, const int* seeds, int count);

#endif // GRAPHENE_UTILS_H

