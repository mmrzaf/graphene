// apps/gwo_im.c
// Grey Wolf Optimization for Influence Maximization (simple, pragmatic)
// Build: make -C apps gwo_im
// Usage example:
//   ./gwo_im ../examples/small.edgelist 5 0.01 50 100 123 50
//   -> <graph> <k> <p> <pop_n> <max_iter> <seed> <ic_runs>

#include "graph.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// ----------------- helpers -----------------

static void fatal(const char *msg) {
  fprintf(stderr, "ERROR: %s\n", msg);
  exit(1);
}

static int cmp_desc_double_idx(const void *a, const void *b) {
  const double *da = (const double *)a;
  const double *db = (const double *)b;
  if (da[0] < db[0])
    return 1;
  if (da[0] > db[0])
    return -1;
  // tie-break by index (stored in second slot)
  if ((int)da[1] < (int)db[1])
    return -1;
  if ((int)da[1] > (int)db[1])
    return 1;
  return 0;
}

// pick top-k candidate indices from position vector pos (size m)
// returns newly allocated array of size k (node indices are indices into
// candidates[])
static int *top_k_from_pos(const double *pos, int m, int k,
                           const int *candidates) {
  if (k <= 0)
    return NULL;
  if (k > m)
    k = m;

  // build array of (value, idx)
  double *arr = (double *)malloc(sizeof(double) * 2 * m);
  if (!arr)
    fatal("malloc failed in top_k_from_pos");
  for (int i = 0; i < m; ++i) {
    arr[2 * i + 0] = pos[i];
    arr[2 * i + 1] = (double)i; // store index as double for sorting tie-break
  }

  // qsort descending
  qsort(arr, m, 2 * sizeof(double), cmp_desc_double_idx);

  int *seeds = (int *)malloc(sizeof(int) * k);
  if (!seeds)
    fatal("malloc failed in top_k_from_pos seeds");
  for (int i = 0; i < k; ++i) {
    int cand_index = (int)arr[2 * i + 1];
    seeds[i] = candidates[cand_index];
  }

  free(arr);
  return seeds;
}

// compute IC spread (Monte Carlo average)
// graph g, seed array (node ids), k seeds, prob p, runs monte_carlo_runs
static double ic_spread(const Graph *g, const int *seeds, int k, double p,
                        int runs) {
  if (!g)
    return 0.0;
  int n = g->n;
  long long total_activated = 0;

  for (int run = 0; run < runs; ++run) {
    // active flags
    char *active = (char *)calloc(n, 1);
    if (!active)
      fatal("malloc failed in ic_spread active");
    int *queue = (int *)malloc(sizeof(int) * n);
    if (!queue)
      fatal("malloc failed in ic_spread queue");
    int qh = 0, qt = 0;

    for (int i = 0; i < k; ++i) {
      int v = seeds[i];
      if (v < 0 || v >= n)
        continue;
      if (!active[v]) {
        active[v] = 1;
        queue[qt++] = v;
      }
    }

    while (qh < qt) {
      int u = queue[qh++];
      int deg;
      const int *nbrs = graph_neighbors(g, u, &deg);
      for (int i = 0; i < deg; ++i) {
        int v = nbrs[i];
        if (v < 0 || v >= n)
          continue;
        if (active[v])
          continue;
        double r = utils_rand_double();
        if (r < p) {
          active[v] = 1;
          queue[qt++] = v;
        }
      }
    }

    int activated = 0;
    for (int i = 0; i < n; ++i)
      activated += active[i];
    total_activated += activated;

    free(active);
    free(queue);
  }

  return (double)total_activated / (double)runs;
}

// compute worthiness and fitness for set S (node ids), using p and graph
// Returns fitness value (higher better). Also returns sum_w via out_sum_w (if
// non-null). Implementation notes:
//  - Prob_receive approximated by combining direct contributions (p) and
//  two-hop (p^2) per path.
//  - Worthiness W(u) = P_receive(u) * deg(u).
//  - Fitness = (normalized Shannon entropy over R(S)) * (sum_w)
//    This encourages sets that reach many worthy nodes and do so evenly.
static double compute_fitness(const Graph *g, const int *S, int k, double p,
                              double *out_sum_w) {
  if (!g)
    return 0.0;
  int n = g->n;

  // mark seeds
  char *is_seed = (char *)calloc(n, 1);
  if (!is_seed)
    fatal("malloc failed in compute_fitness is_seed");
  for (int i = 0; i < k; ++i) {
    if (S[i] >= 0 && S[i] < n)
      is_seed[S[i]] = 1;
  }

  // We'll accumulate contributions per node as product complement:
  // prob_receive = 1 - prod(1 - contrib) where contrib is p for direct edge,
  // p^2 for 2-hop path
  double *one_minus = (double *)malloc(sizeof(double) * n);
  if (!one_minus)
    fatal("malloc failed in compute_fitness one_minus");
  for (int i = 0; i < n; ++i)
    one_minus[i] = 1.0;

  // For all seeds, add direct neighbors contributions and second-hop via
  // neighbors
  for (int i = 0; i < k; ++i) {
    int s = S[i];
    if (s < 0 || s >= n)
      continue;
    int deg_s;
    const int *nbrs_s = graph_neighbors(g, s, &deg_s);
    // direct neighbors (contrib p)
    for (int j = 0; j < deg_s; ++j) {
      int u = nbrs_s[j];
      // if seed itself, skip (we can mark it later)
      double contrib = p;
      one_minus[u] *= (1.0 - contrib);
      // second hop: neighbors of u get contrib p^2 (approx)
      int deg_u;
      const int *nbrs_u = graph_neighbors(g, u, &deg_u);
      for (int t = 0; t < deg_u; ++t) {
        int v = nbrs_u[t];
        // avoid adding for the seed node s or immediate duplicates
        if (v == s)
          continue;
        double contrib2 = p * p;
        one_minus[v] *= (1.0 - contrib2);
      }
    }
    // also consider seed itself receives probability 1 (by definition)
    one_minus[s] *= 0.0; // set to zero => p_receive = 1.0
  }

  // collect R(S): nodes with any change (one_minus != 1) or neighbors of seeds
  char *in_R = (char *)calloc(n, 1);
  if (!in_R)
    fatal("malloc failed in compute_fitness in_R");
  int Rcount = 0;
  for (int i = 0; i < k; ++i) {
    int s = S[i];
    if (s < 0 || s >= n)
      continue;
    if (!in_R[s]) {
      in_R[s] = 1;
      Rcount++;
    }
    int deg_s;
    const int *nbrs_s = graph_neighbors(g, s, &deg_s);
    for (int j = 0; j < deg_s; ++j) {
      int u = nbrs_s[j];
      if (!in_R[u]) {
        in_R[u] = 1;
        Rcount++;
      }
      // second neighbors
      int deg_u;
      const int *nbrs_u = graph_neighbors(g, u, &deg_u);
      for (int t = 0; t < deg_u; ++t) {
        int v = nbrs_u[t];
        if (!in_R[v]) {
          in_R[v] = 1;
          Rcount++;
        }
      }
    }
  }

  if (Rcount == 0) {
    free(is_seed);
    free(one_minus);
    free(in_R);
    if (out_sum_w)
      *out_sum_w = 0.0;
    return 0.0;
  }

  // compute worthiness W(u) = p_receive * deg(u)
  double *W = (double *)malloc(sizeof(double) * Rcount);
  if (!W)
    fatal("malloc failed in compute_fitness W");
  int idx = 0;
  double sumW = 0.0;
  for (int u = 0; u < n; ++u) {
    if (!in_R[u])
      continue;
    double p_receive = 1.0 - one_minus[u];
    int deg_u = graph_degree(g, u);
    double w = p_receive * (double)deg_u;
    W[idx++] = w;
    sumW += w;
  }

  // If sumW == 0, assign tiny positive to avoid log(0)
  if (sumW <= 0.0) {
    for (int i = 0; i < Rcount; ++i)
      W[i] = 1e-12;
    sumW = 1e-12 * Rcount;
  }

  // compute normalized Shannon entropy H = -sum p_i ln p_i
  double H = 0.0;
  for (int i = 0; i < Rcount; ++i) {
    double pi = W[i] / sumW;
    if (pi > 0.0)
      H -= pi * log(pi);
  }
  // Normalize entropy by log(Rcount) to get [0,1]
  double Hnorm = (Rcount > 1) ? (H / log((double)Rcount)) : 1.0;

  // final fitness: Hnorm * sumW (simple mixing of evenness and magnitude)
  double fitness = Hnorm * sumW;

  if (out_sum_w)
    *out_sum_w = sumW;

  free(is_seed);
  free(one_minus);
  free(in_R);
  free(W);
  return fitness;
}

// map: build candidate list (exclude degree 1 nodes)
static int *build_candidates(const Graph *g, int *out_m) {
  int n = g->n;
  int *tmp = (int *)malloc(sizeof(int) * n);
  if (!tmp)
    fatal("malloc failed in build_candidates tmp");
  int m = 0;
  for (int v = 0; v < n; ++v) {
    int deg = graph_degree(g, v);
    if (deg > 1) {
      tmp[m++] = v;
    }
  }
  // shrink
  int *cands = (int *)malloc(sizeof(int) * m);
  if (!cands)
    fatal("malloc failed in build_candidates cands");
  memcpy(cands, tmp, m * sizeof(int));
  free(tmp);
  *out_m = m;
  return cands;
}

// ----------------- GWO algorithm -----------------

int main(int argc, char **argv) {
  if (argc < 8) {
    fprintf(stderr,
            "Usage: %s <graph.edgelist> <k> <p> <pop_n> <max_iter> <seed> "
            "<ic_runs>\n"
            "Example: %s ../examples/small.edgelist 5 0.01 50 100 123 50\n",
            argv[0], argv[0]);
    return 1;
  }

  const char *path = argv[1];
  int k = atoi(argv[2]);
  double p = atof(argv[3]);
  int pop_n = atoi(argv[4]);
  int max_iter = atoi(argv[5]);
  unsigned int seed = (unsigned int)atoi(argv[6]);
  int ic_runs = atoi(argv[7]);

  if (k <= 0)
    fatal("k must be > 0");
  if (p <= 0.0 || p > 1.0)
    fatal("p must be in (0,1]");
  if (pop_n <= 0)
    pop_n = 50;
  if (max_iter <= 0)
    max_iter = 100;
  if (ic_runs <= 0)
    ic_runs = 50;

  utils_seed_rng(seed);

  // load graph (keep it as undirected for influence, if possible)
  Graph *g = graph_load_with_flags(path, GRAPH_FLAG_NONE);
  if (!g) {
    fprintf(stderr, "graph_load failed: %s\n", graph_last_error());
    return 1;
  }

  int n = g->n;
  printf("Loaded graph: n=%d m=%d\n", n, graph_num_edges(g));

  // build candidate pool (degree > 1)
  int m;
  int *candidates = build_candidates(g, &m);
  printf("Candidate pool size (deg>1) m=%d\n", m);
  if (m == 0)
    fatal("No eligible candidates (all nodes degree <= 1)");

  if (k > m) {
    fprintf(stderr,
            "Warning: k (%d) > candidate count (%d). Reducing k to %d\n", k, m,
            m);
    k = m;
  }

  // degrees of candidates (for seeding initialization)
  double *deg_c = (double *)malloc(sizeof(double) * m);
  if (!deg_c)
    fatal("malloc failed deg_c");
  for (int i = 0; i < m; ++i)
    deg_c[i] = (double)graph_degree(g, candidates[i]);

  // allocate population: positions (pop_n x m)
  double **X = (double **)malloc(sizeof(double *) * pop_n);
  double *fitness = (double *)malloc(sizeof(double) * pop_n);
  double *sumW = (double *)malloc(sizeof(double) * pop_n);
  if (!X || !fitness || !sumW)
    fatal("malloc failed population");

  for (int i = 0; i < pop_n; ++i) {
    X[i] = (double *)malloc(sizeof(double) * m);
    if (!X[i])
      fatal("malloc failed X[i]");
  }

  // init wolves: bias by degree * rand
  for (int i = 0; i < pop_n; ++i) {
    for (int j = 0; j < m; ++j) {
      X[i][j] = deg_c[j] * (0.5 + utils_rand_double()); // bias but keep spread
      if (X[i][j] < 0.0)
        X[i][j] = 0.0;
    }
    // compute fitness via top-k
    int *seeds = top_k_from_pos(X[i], m, k, candidates);
    fitness[i] = compute_fitness(g, seeds, k, p, &sumW[i]);
    free(seeds);
  }

  // identify alpha/beta/delta
  int alpha = 0, beta = 0, delta = 0;
  for (int it = 0; it < max_iter; ++it) {
    // find best three
    alpha = beta = delta = 0;
    for (int i = 0; i < pop_n; ++i) {
      if (fitness[i] > fitness[alpha])
        alpha = i;
    }
    for (int i = 0; i < pop_n; ++i) {
      if (i == alpha)
        continue;
      if (fitness[i] > fitness[beta])
        beta = i;
    }
    for (int i = 0; i < pop_n; ++i) {
      if (i == alpha || i == beta)
        continue;
      if (fitness[i] > fitness[delta])
        delta = i;
    }

    // linearly decrease a from 2 -> 0
    double a_param = 2.0 * (1.0 - ((double)it / (double)MAX(1, max_iter - 1)));

    // Update positions for omegas (all wolves)
    for (int w = 0; w < pop_n; ++w) {
      if (w == alpha)
        continue; // keep alpha
      // compute new Xw
      for (int j = 0; j < m; ++j) {
        double r1 = utils_rand_double();
        double r2 = utils_rand_double();
        double A1 = 2.0 * a_param * r1 - a_param;
        double C1 = 2.0 * r2;

        double D_alpha = fabs(C1 * X[alpha][j] - X[w][j]);
        double X1 = X[alpha][j] - A1 * D_alpha;

        r1 = utils_rand_double();
        r2 = utils_rand_double();
        double A2 = 2.0 * a_param * r1 - a_param;
        double C2 = 2.0 * r2;
        double D_beta = fabs(C2 * X[beta][j] - X[w][j]);
        double X2 = X[beta][j] - A2 * D_beta;

        r1 = utils_rand_double();
        r2 = utils_rand_double();
        double A3 = 2.0 * a_param * r1 - a_param;
        double C3 = 2.0 * r2;
        double D_delta = fabs(C3 * X[delta][j] - X[w][j]);
        double X3 = X[delta][j] - A3 * D_delta;

        double newx = (X1 + X2 + X3) / 3.0;
        if (newx < 0.0)
          newx = 0.0; // clip to non-negative
        X[w][j] = newx;
      }

      // rebuild seeds and evaluate fitness
      int *seeds = top_k_from_pos(X[w], m, k, candidates);
      double new_fit = compute_fitness(g, seeds, k, p, &sumW[w]);
      free(seeds);
      fitness[w] = new_fit;
    }

    // optional regeneration: if beta or delta identical to alpha (by fitness),
    // randomize them
    if (fitness[beta] == fitness[alpha]) {
      for (int j = 0; j < m; ++j)
        X[beta][j] = deg_c[j] * (0.5 + utils_rand_double());
      int *s = top_k_from_pos(X[beta], m, k, candidates);
      fitness[beta] = compute_fitness(g, s, k, p, &sumW[beta]);
      free(s);
    }
    if (fitness[delta] == fitness[alpha] || fitness[delta] == fitness[beta]) {
      for (int j = 0; j < m; ++j)
        X[delta][j] = deg_c[j] * (0.5 + utils_rand_double());
      int *s = top_k_from_pos(X[delta], m, k, candidates);
      fitness[delta] = compute_fitness(g, s, k, p, &sumW[delta]);
      free(s);
    }

    if (it % 10 == 0 || it == max_iter - 1) {
      printf("iter %d: best fitness=%.6f (alpha=%d)\n", it, fitness[alpha],
             alpha);
    }
  }

  // after iterations, pick alpha's seed set
  int *best_seeds = top_k_from_pos(X[alpha], m, k, candidates);
  printf("Best seed set (alpha): ");
  for (int i = 0; i < k; ++i)
    printf("%d ", best_seeds[i]);
  printf("\n");

  // final evaluation with IC Monte Carlo
  double avg_spread = ic_spread(g, best_seeds, k, p, ic_runs);
  printf("Final IC spread (avg over %d runs): %.4f nodes\n", ic_runs,
         avg_spread);

  // cleanup
  free(best_seeds);
  for (int i = 0; i < pop_n; ++i)
    free(X[i]);
  free(X);
  free(fitness);
  free(sumW);
  free(deg_c);
  free(candidates);
  graph_free(g);
  return 0;
}
