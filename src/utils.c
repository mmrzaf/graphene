#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static unsigned int rng_state = 1;

void utils_seed_rng(unsigned int seed) { rng_state = seed; }

int utils_rand_int(int min, int max) {
  if (min > max) {
    int temp = min;
    min = max;
    max = temp;
  }

  rng_state = rng_state * 1103515245 + 12345;
  unsigned int val = (rng_state / 65536) % 32768;
  return min + (val % (max - min + 1));
}

double utils_rand_double(void) {
  rng_state = rng_state * 1103515245 + 12345;
  unsigned int val = (rng_state / 65536) % 32768;
  return val / 32768.0;
}

void utils_shuffle_int_array(int *arr, int n) {
  for (int i = n - 1; i > 0; i--) {
    int j = utils_rand_int(0, i);
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
  }
}

int *utils_load_seedset(const char *path, int *out_count) {
  FILE *f = fopen(path, "r");
  if (!f) {
    if (out_count)
      *out_count = 0;
    return NULL;
  }

  int capacity = 16;
  int count = 0;
  int *seeds = (int *)malloc(capacity * sizeof(int));

  int node;
  while (fscanf(f, "%d", &node) == 1) {
    if (count >= capacity) {
      capacity *= 2;
      int *new_seeds = (int *)realloc(seeds, capacity * sizeof(int));
      if (!new_seeds) {
        free(seeds);
        fclose(f);
        if (out_count)
          *out_count = 0;
        return NULL;
      }
      seeds = new_seeds;
    }
    seeds[count++] = node;
  }

  fclose(f);

  if (out_count)
    *out_count = count;
  return seeds;
}

void utils_save_seedset(const char *path, const int *seeds, int count) {
  FILE *f = fopen(path, "w");
  if (!f)
    return;

  for (int i = 0; i < count; i++) {
    fprintf(f, "%d\n", seeds[i]);
  }

  fclose(f);
}
