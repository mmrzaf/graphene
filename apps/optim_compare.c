// apps/optim_compare.c
// Compare 3 local-improvement strategies in GWO for Influence Maximization.
// Build: make -C apps optim_compare
// Usage:
//  ./optim_compare <graph.edgelist> <k> <p> <pop_n> <max_iter> <seed> <ic_runs>
//  <out_prefix>
// Example:
//  ./optim_compare ../examples/small.edgelist 5 0.01 40 100 123 200
//  results/run1

#define _POSIX_C_SOURCE 200809L
#include "graph.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

// ------------ Helpers (adapted and simplified) ------------

static void fatal(const char *msg) {
  fprintf(stderr, "FATAL: %s\n", msg);
  exit(1);
}

static int int_cmp_inc(const void *a, const void *b) {
  int x = *(const int *)a;
  int y = *(const int *)b;
  return (x > y) - (x < y);
}

// top-k selectors: pos is size m (scores for candidates[]), returns allocated
// array length k
static int *top_k_from_pos(const double *pos, int m, int k,
                           const int *candidates) {
  if (k <= 0)
    return NULL;
  if (k > m)
    k = m;
  // array of indices 0..m-1
  int *idx = malloc(sizeof(int) * m);
  if (!idx)
    fatal("malloc");
  for (int i = 0; i < m; ++i)
    idx[i] = i;
  // partial selection using simple sort since m likely moderate; easier to
  // reason about
  qsort(idx, m, sizeof(int), (int (*)(const void *, const void *))int_cmp_inc);
  // but idx is ascending; we want top descending, so pick last k
  int *seeds = malloc(sizeof(int) * k);
  if (!seeds)
    fatal("malloc");
  for (int i = 0; i < k; ++i) {
    int pick = idx[m - 1 - i];
    seeds[i] = candidates[pick];
  }
  free(idx);
  return seeds;
}

// Monte Carlo IC spread (average nodes activated)
static double ic_spread(const Graph *g, const int *seeds, int k, double p,
                        int runs) {
  if (!g)
    return 0.0;
  int n = g->n;
  long long total_activated = 0;
  for (int run = 0; run < runs; ++run) {
    char *active = calloc(n, 1);
    int *queue = malloc(sizeof(int) * n);
    if (!active || !queue)
      fatal("malloc");
    int qh = 0, qt = 0;
    for (int i = 0; i < k; ++i) {
      int v = seeds[i];
      if (v >= 0 && v < n && !active[v]) {
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
        if (utils_rand_double() < p) {
          active[v] = 1;
          queue[qt++] = v;
        }
      }
    }
    int cnt = 0;
    for (int i = 0; i < n; ++i)
      cnt += active[i];
    total_activated += cnt;
    free(active);
    free(queue);
  }
  return (double)total_activated / (double)runs;
}

// Simple approximate fitness function (same spirit as earlier): worthiness from
// 1- and 2-hop contributions
static double compute_fitness(const Graph *g, const int *S, int k, double p,
                              double *out_sumw) {
  if (!g)
    return 0.0;
  int n = g->n;
  char *is_seed = calloc(n, 1);
  if (!is_seed)
    fatal("malloc");
  for (int i = 0; i < k; ++i)
    if (S[i] >= 0 && S[i] < n)
      is_seed[S[i]] = 1;

  double *one_minus = malloc(sizeof(double) * n);
  if (!one_minus)
    fatal("malloc");
  for (int i = 0; i < n; ++i)
    one_minus[i] = 1.0;

  // acc contributions
  for (int i = 0; i < k; ++i) {
    int s = S[i];
    if (s < 0 || s >= n)
      continue;
    int degs;
    const int *nbrs = graph_neighbors(g, s, &degs);
    for (int j = 0; j < degs; ++j) {
      int u = nbrs[j];
      one_minus[u] *= (1.0 - p);
      // second hop
      int degu;
      const int *nbrsu = graph_neighbors(g, u, &degu);
      for (int t = 0; t < degu; ++t) {
        int v = nbrsu[t];
        if (v == s)
          continue;
        one_minus[v] *= (1.0 - p * p);
      }
    }
    // seed itself: set receive prob to 1
    one_minus[s] = 0.0;
  }

  // collect R set
  char *inR = calloc(n, 1);
  if (!inR)
    fatal("malloc");
  int Rcount = 0;
  for (int i = 0; i < k; ++i) {
    int s = S[i];
    if (s < 0 || s >= n)
      continue;
    if (!inR[s]) {
      inR[s] = 1;
      Rcount++;
    }
    int degs;
    const int *nbrs = graph_neighbors(g, s, &degs);
    for (int j = 0; j < degs; ++j) {
      int u = nbrs[j];
      if (!inR[u]) {
        inR[u] = 1;
        Rcount++;
      }
      int degu;
      const int *nbrsu = graph_neighbors(g, u, &degu);
      for (int t = 0; t < degu; ++t) {
        int v = nbrsu[t];
        if (!inR[v]) {
          inR[v] = 1;
          Rcount++;
        }
      }
    }
  }

  if (Rcount == 0) {
    free(is_seed);
    free(one_minus);
    free(inR);
    if (out_sumw)
      *out_sumw = 0.0;
    return 0.0;
  }

  double *W = malloc(sizeof(double) * Rcount);
  if (!W)
    fatal("malloc");
  int idx = 0;
  double sumW = 0.0;
  for (int u = 0; u < n; ++u) {
    if (!inR[u])
      continue;
    double precv = 1.0 - one_minus[u];
    int degu = graph_degree(g, u);
    double w = precv * (double)degu;
    W[idx++] = w;
    sumW += w;
  }

  if (sumW <= 0.0) {
    for (int i = 0; i < Rcount; ++i)
      W[i] = 1e-12;
    sumW = 1e-12 * Rcount;
  }

  double H = 0.0;
  for (int i = 0; i < Rcount; ++i) {
    double pi = W[i] / sumW;
    if (pi > 0.0)
      H -= pi * log(pi);
  }
  double Hnorm = (Rcount > 1) ? (H / log((double)Rcount)) : 1.0;
  double fitness = Hnorm * sumW;

  if (out_sumw)
    *out_sumw = sumW;
  free(is_seed);
  free(one_minus);
  free(inR);
  free(W);
  return fitness;
}

// build candidate pool (exclude deg <= 1)
static int *build_candidates(const Graph *g, int *out_m) {
  int n = g->n;
  int *tmp = malloc(sizeof(int) * n);
  if (!tmp)
    fatal("malloc");
  int m = 0;
  for (int v = 0; v < n; ++v) {
    int d = graph_degree(g, v);
    if (d > 1)
      tmp[m++] = v;
  }
  int *cands = malloc(sizeof(int) * m);
  if (!cands && m > 0)
    fatal("malloc");
  memcpy(cands, tmp, m * sizeof(int));
  free(tmp);
  *out_m = m;
  return cands;
}

// convert seed array to JSON array string (caller frees returned char*)
static char *seedset_to_json(const int *S, int k) {
  size_t cap = 64 + k * 8;
  char *buf = malloc(cap);
  if (!buf)
    fatal("malloc");
  size_t pos = 0;
  pos += snprintf(buf + pos, cap - pos, "[");
  for (int i = 0; i < k; ++i) {
    pos += snprintf(buf + pos, cap - pos, "%d", S[i]);
    if (i + 1 < k)
      pos += snprintf(buf + pos, cap - pos, ",");
  }
  pos += snprintf(buf + pos, cap - pos, "]");
  return buf;
}

// ------------ Local optimization rules ------------

// Try neighbor replacement: for each seed s, examine neighbors of s and pick
// neighbor candidate that improves fitness (test neighbors sequentially; accept
// first improving move). Records action via action_buf (allocated string),
// caller must free. Returns 1 if any change applied, 0 otherwise.
static int local_optimize_neighbors(const Graph *g, int *S, int k, double p,
                                    const int *candidates, int m,
                                    double *out_newfit, char **action_buf) {
  double base_sumw = 0.0;
  double base_fit = compute_fitness(g, S, k, p, &base_sumw);
  int changed = 0;
  // for fast membership test
  int n = g->n;
  char *ismask = calloc(n, 1);
  if (!ismask)
    fatal("malloc");
  for (int i = 0; i < k; ++i)
    if (S[i] >= 0 && S[i] < n)
      ismask[S[i]] = 1;

  // buffer for action description
  char *buf = malloc(1024);
  if (!buf)
    fatal("malloc");
  buf[0] = '\0';
  size_t bpos = 0;
  bpos += snprintf(buf + bpos, 1024 - bpos,
                   "{\"type\":\"neighbors\",\"base_fit\":%.6f,\"changes\":[",
                   base_fit);

  for (int i = 0; i < k; ++i) {
    int s = S[i];
    if (s < 0 || s >= n)
      continue;
    int deg;
    const int *nbrs = graph_neighbors(g, s, &deg);
    for (int j = 0; j < deg; ++j) {
      int cand = nbrs[j];
      if (cand < 0 || cand >= n)
        continue;
      if (ismask[cand])
        continue; // already in seed set
      // propose replacing S[i] with cand
      int old = S[i];
      S[i] = cand;
      double new_sumw = 0.0;
      double newfit = compute_fitness(g, S, k, p, &new_sumw);
      if (newfit > base_fit + 1e-12) {
        // accept
        changed = 1;
        base_fit = newfit;
        bpos += snprintf(buf + bpos, 1024 - bpos,
                         "{\"slot\":%d,\"from\":%d,\"to\":%d,\"newfit\":%.6f},",
                         i, old, cand, newfit);
        // update mask
        ismask[old] = 0;
        ismask[cand] = 1;
        break; // move to next seed
      } else {
        // revert
        S[i] = old;
      }
    }
  }

  // trim trailing comma if present and close JSON
  if (bpos > 0 && buf[bpos - 1] == ',')
    buf[bpos - 1] = '\0';
  bpos += snprintf(buf + bpos, 1024 - bpos, "]}");
  *action_buf = buf;
  free(ismask);
  *out_newfit = base_fit;
  return changed;
}

// BFS shortest path from s to t, returns path array (allocated) and length via
// out_len. If no path, returns NULL and out_len=0.
static int *bfs_path(const Graph *g, int s, int t, int *out_len) {
  int n = g->n;
  char *vis = calloc(n, 1);
  int *parent = malloc(sizeof(int) * n);
  int *queue = malloc(sizeof(int) * n);
  if (!vis || !parent || !queue)
    fatal("malloc");
  for (int i = 0; i < n; i++)
    parent[i] = -1;
  int qh = 0, qt = 0;
  queue[qt++] = s;
  vis[s] = 1;
  parent[s] = -1;
  int found = 0;
  while (qh < qt) {
    int u = queue[qh++];
    if (u == t) {
      found = 1;
      break;
    }
    int deg;
    const int *nbrs = graph_neighbors(g, u, &deg);
    for (int i = 0; i < deg; i++) {
      int v = nbrs[i];
      if (!vis[v]) {
        vis[v] = 1;
        parent[v] = u;
        queue[qt++] = v;
      }
    }
  }
  if (!found) {
    free(vis);
    free(parent);
    free(queue);
    *out_len = 0;
    return NULL;
  }
  // reconstruct
  int cur = t;
  int pathlen = 0;
  while (cur != -1) {
    pathlen++;
    cur = parent[cur];
  }
  int *path = malloc(sizeof(int) * pathlen);
  cur = t;
  for (int i = pathlen - 1; i >= 0; --i) {
    path[i] = cur;
    cur = parent[cur];
  }
  free(vis);
  free(parent);
  free(queue);
  *out_len = pathlen;
  return path;
}

// For each pair of seeds, consider internal nodes on shortest path between them
// (excluding endpoints). Try replacing any seed by these internal nodes if
// fitness improves. Record actions in JSON.
static int local_optimize_shortestpath(const Graph *g, int *S, int k, double p,
                                       const int *candidates, int m,
                                       double *out_newfit, char **action_buf) {
  double base_sumw = 0.0;
  double base_fit = compute_fitness(g, S, k, p, &base_sumw);
  int n = g->n;
  char *ismask = calloc(n, 1);
  if (!ismask)
    fatal("malloc");
  for (int i = 0; i < k; ++i)
    if (S[i] >= 0 && S[i] < n)
      ismask[S[i]] = 1;

  char *buf = malloc(2048);
  if (!buf)
    fatal("malloc");
  buf[0] = '\0';
  size_t bpos = 0;
  bpos += snprintf(buf + bpos, 2048 - bpos,
                   "{\"type\":\"shortestpath\",\"base_fit\":%.6f,\"changes\":[",
                   base_fit);

  int changed = 0;
  for (int a = 0; a < k; ++a) {
    for (int b = a + 1; b < k; ++b) {
      int s = S[a], t = S[b];
      if (s < 0 || t < 0 || s >= n || t >= n)
        continue;
      int plen = 0;
      int *path = bfs_path(g, s, t, &plen);
      if (!path || plen <= 2) {
        free(path);
        continue;
      } // no internal nodes
      for (int pi = 1; pi < plen - 1; ++pi) {
        int cand = path[pi];
        if (ismask[cand])
          continue;
        // try replacing seed a with cand
        int old = S[a];
        S[a] = cand;
        double new_sum = 0.0;
        double newfit = compute_fitness(g, S, k, p, &new_sum);
        if (newfit > base_fit + 1e-12) {
          changed = 1;
          base_fit = newfit;
          bpos +=
              snprintf(buf + bpos, 2048 - bpos,
                       "{\"slot\":%d,\"from\":%d,\"to\":%d,\"newfit\":%.6f},",
                       a, old, cand, newfit);
          ismask[old] = 0;
          ismask[cand] = 1;
          break;
        } else {
          S[a] = old;
        }
      }
      free(path);
    }
  }

  if (bpos > 0 && buf[bpos - 1] == ',')
    buf[bpos - 1] = '\0';
  bpos += snprintf(buf + bpos, 2048 - bpos, "]}");
  *action_buf = buf;
  free(ismask);
  *out_newfit = base_fit;
  return changed;
}

// ------------ GWO main runner (single-mode) ------------

// mode: 1 = no optimization; 2 = neighbor replacement; 3 = shortest-path
// replacement
static void run_mode(const Graph *g, int k, double p, int pop_n, int max_iter,
                     int seed, int ic_runs, int mode, const int *candidates,
                     int m, const char *out_prefix) {
  char logfile[1024];
  snprintf(logfile, sizeof(logfile), "%s_mode%d.jsonl", out_prefix, mode);
  FILE *logf = fopen(logfile, "w");
  if (!logf) {
    fprintf(
        stderr,
        "Warning: cannot open log file %s, proceeding without file logging.\n",
        logfile);
  }

  utils_seed_rng(seed + mode); // slight variation per mode

  int n = g->n;
  double **X = malloc(sizeof(double *) * pop_n);
  double *fit = malloc(sizeof(double) * pop_n);
  double *sumw = malloc(sizeof(double) * pop_n);
  if (!X || !fit || !sumw)
    fatal("malloc");
  for (int i = 0; i < pop_n; ++i) {
    X[i] = malloc(sizeof(double) * m);
    if (!X[i])
      fatal("malloc");
  }

  // initialize by degree-bias
  double *degc = malloc(sizeof(double) * m);
  if (!degc)
    fatal("malloc");
  for (int i = 0; i < m; ++i)
    degc[i] = (double)graph_degree(g, candidates[i]);

  for (int i = 0; i < pop_n; ++i) {
    for (int j = 0; j < m; ++j)
      X[i][j] = degc[j] * (0.5 + utils_rand_double());
    int *s = top_k_from_pos(X[i], m, k, candidates);
    fit[i] = compute_fitness(g, s, k, p, &sumw[i]);
    free(s);
  }

  int alpha = 0, beta = 0, delta = 0;
  for (int it = 0; it < max_iter; ++it) {
    // identify top3
    alpha = beta = delta = 0;
    for (int i = 0; i < pop_n; ++i)
      if (fit[i] > fit[alpha])
        alpha = i;
    for (int i = 0; i < pop_n; ++i)
      if (i != alpha && fit[i] > fit[beta])
        beta = i;
    for (int i = 0; i < pop_n; ++i)
      if (i != alpha && i != beta && fit[i] > fit[delta])
        delta = i;

    double a_param = 2.0 * (1.0 - ((double)it / (double)MAX(1, max_iter - 1)));

    for (int w = 0; w < pop_n; ++w) {
      if (w == alpha)
        continue;
      for (int j = 0; j < m; ++j) {
        double r1 = utils_rand_double(), r2 = utils_rand_double();
        double A1 = 2.0 * a_param * r1 - a_param;
        double C1 = 2.0 * r2;
        double D1 = fabs(C1 * X[alpha][j] - X[w][j]);
        double X1 = X[alpha][j] - A1 * D1;

        r1 = utils_rand_double();
        r2 = utils_rand_double();
        double A2 = 2.0 * a_param * r1 - a_param;
        double C2 = 2.0 * r2;
        double D2 = fabs(C2 * X[beta][j] - X[w][j]);
        double X2 = X[beta][j] - A2 * D2;

        r1 = utils_rand_double();
        r2 = utils_rand_double();
        double A3 = 2.0 * a_param * r1 - a_param;
        double C3 = 2.0 * r2;
        double D3 = fabs(C3 * X[delta][j] - X[w][j]);
        double X3 = X[delta][j] - A3 * D3;

        double newx = (X1 + X2 + X3) / 3.0;
        if (newx < 0.0)
          newx = 0.0;
        X[w][j] = newx;
      }

      // get seed set for this wolf
      int *seeds = top_k_from_pos(X[w], m, k, candidates);
      double before_fit = compute_fitness(g, seeds, k, p, NULL);
      double after_fit = before_fit;
      char *action = NULL;

      // local optimization according to mode
      if (mode == 2) {
        int changed = local_optimize_neighbors(g, seeds, k, p, candidates, m,
                                               &after_fit, &action);
        (void)changed;
      } else if (mode == 3) {
        int changed = local_optimize_shortestpath(g, seeds, k, p, candidates, m,
                                                  &after_fit, &action);
        (void)changed;
      } else {
        // mode 1 -> no optimization; produce minimal action record
        action = malloc(128);
        snprintf(action, 128, "{\"type\":\"none\"}");
      }

      // if improved, update X[w] positions to reflect new seeds:
      if (after_fit > before_fit + 1e-12) {
        // we 'nudge' X[w] so that chosen candidate indices get slightly higher
        // values Build an index map from candidate node -> pos index
        for (int j = 0; j < m; ++j)
          X[w][j] *= 0.95; // slight shrink
        for (int sidx = 0; sidx < k; ++sidx) {
          int node = seeds[sidx];
          // find index in candidates
          for (int j = 0; j < m; ++j) {
            if (candidates[j] == node) {
              X[w][j] += 1.0;
              break;
            }
          }
        }
      }

      // update fitness record for wolf
      fit[w] = after_fit;

      // log iteration: JSON line
      if (logf) {
        char *before_js = seedset_to_json(seeds, k);
        // need to write new seed set (seeds may have been changed by local opt)
        char *after_js = seedset_to_json(seeds, k);
        time_t now = time(NULL);
        fprintf(logf,
                "{\"time\":%ld,\"iter\":%d,\"wolf\":%d,\"mode\":%d,\"before_"
                "fit\":%.6f,"
                "\"after_fit\":%.6f,\"before_seeds\":%s,\"after_seeds\":%s,"
                "\"action\":%s}\n",
                (long)now, it, w, mode, before_fit, after_fit, before_js,
                after_js, action);
        fflush(logf);
        free(before_js);
        free(after_js);
      }

      free(action);
      free(seeds);
    } // end wolves

    if (it % 10 == 0 || it == max_iter - 1) {
      fprintf(stderr, "[mode %d] iter %d alpha_fit=%.6f\n", mode, it,
              fit[alpha]);
    }
  } // end iterations

  // finalize: pick alpha seeds
  // recompute alpha
  alpha = 0;
  for (int i = 0; i < pop_n; ++i)
    if (fit[i] > fit[alpha])
      alpha = i;
  int *best_seeds = top_k_from_pos(X[alpha], m, k, candidates);
  double final_fit = compute_fitness(g, best_seeds, k, p, NULL);
  double final_spread = ic_spread(g, best_seeds, k, p, ic_runs);

  // summary JSON
  char summfile[1024];
  snprintf(summfile, sizeof(summfile), "%s_mode%d_summary.json", out_prefix,
           mode);
  FILE *sf = fopen(summfile, "w");
  if (sf) {
    char *js = seedset_to_json(best_seeds, k);
    fprintf(sf,
            "{\n  \"mode\": %d,\n  \"k\": %d,\n  \"final_fit\": %.6f,\n  "
            "\"final_spread\": %.6f,\n  \"ic_runs\": %d,\n  \"best_seeds\": "
            "%s\n}\n",
            mode, k, final_fit, final_spread, ic_runs, js);
    fclose(sf);
    free(js);
  } else {
    fprintf(stderr, "Could not write summary file %s\n", summfile);
  }

  // print short summary to stdout
  printf(
      "MODE %d SUMMARY: final_fit=%.6f final_spread(avg nodes) = %.4f seeds: ",
      mode, final_fit, final_spread);
  for (int i = 0; i < k; ++i)
    printf("%d ", best_seeds[i]);
  printf("\n");
  if (logf)
    fclose(logf);

  // cleanup
  for (int i = 0; i < pop_n; ++i)
    free(X[i]);
  free(X);
  free(fit);
  free(sumw);
  free(degc);
  free(best_seeds);
}

// ------------ main: run three modes sequentially ------------

int main(int argc, char **argv) {
  if (argc < 9) {
    fprintf(stderr,
            "Usage: %s <graph.edgelist> <k> <p> <pop_n> <max_iter> <seed> "
            "<ic_runs> <out_prefix>\n"
            "Example: %s ../examples/small.edgelist 5 0.01 40 100 123 200 "
            "results/run\n",
            argv[0], argv[0]);
    return 1;
  }

  const char *path = argv[1];
  int k = atoi(argv[2]);
  double p = atof(argv[3]);
  int pop_n = atoi(argv[4]);
  int max_iter = atoi(argv[5]);
  int seed = atoi(argv[6]);
  int ic_runs = atoi(argv[7]);
  const char *out_prefix = argv[8];

  if (k <= 0)
    fatal("k must be > 0");
  if (p <= 0.0 || p > 1.0)
    fatal("p must be (0,1]");
  if (pop_n <= 4)
    pop_n = 20;
  if (max_iter <= 0)
    max_iter = 50;
  if (ic_runs <= 0)
    ic_runs = 100;

  Graph *g = graph_load_with_flags(path, GRAPH_FLAG_NONE);
  if (!g) {
    fprintf(stderr, "graph_load failed: %s\n", graph_last_error());
    return 1;
  }
  printf("Loaded graph: n=%d m=%d\n", g->n, graph_num_edges(g));

  int m;
  int *candidates = build_candidates(g, &m);
  if (m == 0)
    fatal("No candidates (deg<=1 everywhere)");
  if (k > m) {
    fprintf(stderr, "k > candidate count, reducing k to %d\n", m);
    k = m;
  }

  // Run three modes: 1=no opt, 2=neighbor, 3=shortest-path
  for (int mode = 1; mode <= 3; ++mode) {
    printf("\n--- Running mode %d ---\n", mode);
    run_mode(g, k, p, pop_n, max_iter, seed, ic_runs, mode, candidates, m,
             out_prefix);
  }

  free(candidates);
  graph_free(g);
  return 0;
}
