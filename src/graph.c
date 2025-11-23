#include "graph.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char g_last_error[256] = {0};

static void set_error(const char *msg) {
  snprintf(g_last_error, sizeof(g_last_error), "%s", msg);
}

const char *graph_last_error(void) { return g_last_error; }

Graph *graph_create(int num_nodes) {
  if (num_nodes < 0) {
    set_error("Invalid number of nodes");
    return NULL;
  }

  Graph *g = (Graph *)malloc(sizeof(Graph));
  if (!g) {
    set_error("Failed to allocate graph");
    return NULL;
  }

  g->n = num_nodes;
  g->m = 0;
  g->flags = GRAPH_FLAG_NONE;
  g->meta_free = NULL;

  g->offsets = (int *)calloc(num_nodes + 1, sizeof(int));
  g->nodes = (Node *)calloc(num_nodes, sizeof(Node));
  g->edges = NULL;

  if (!g->offsets || !g->nodes) {
    set_error("Allocation failure during graph_create");
    graph_free(g);
    return NULL;
  }

  for (int i = 0; i < num_nodes; i++) {
    g->nodes[i].id = i;
    g->nodes[i].metadata = NULL;
  }

  return g;
}

static int int_cmp(const void *a, const void *b) {
  int x = *(const int *)a;
  int y = *(const int *)b;
  return (x > y) - (x < y);
}

Graph *graph_load(const char *path) {
  return graph_load_with_flags(path, GRAPH_FLAG_NONE);
}

Graph *graph_load_with_flags(const char *path, unsigned int flags) {
  FILE *f = fopen(path, "r");
  if (!f) {
    set_error("Failed to open file");
    return NULL;
  }

  long long u_ll, v_ll;
  long long max_node = -1;
  size_t edge_count = 0;
  int scanned;

  while ((scanned = fscanf(f, "%lld %lld", &u_ll, &v_ll)) != EOF) {
    if (scanned != 2) {
      set_error("Malformed line in input");
      fclose(f);
      return NULL;
    }
    if (u_ll < 0 || v_ll < 0) {
      set_error("Negative node id");
      fclose(f);
      return NULL;
    }
    if (u_ll > max_node)
      max_node = u_ll;
    if (v_ll > max_node)
      max_node = v_ll;
    edge_count++;
    if (flags & GRAPH_FLAG_UNDIRECTED)
      edge_count++;
    if (edge_count > (size_t)INT_MAX) {
      set_error("Edge count overflow");
      fclose(f);
      return NULL;
    }
  }

  int num_nodes = (int)(max_node + 1);
  Graph *g = graph_create(num_nodes);
  if (!g) {
    fclose(f);
    return NULL;
  }
  g->flags = flags;

  size_t *degrees = (size_t *)calloc(num_nodes, sizeof(size_t));
  if (!degrees) {
    set_error("Allocation failure (degrees)");
    graph_free(g);
    fclose(f);
    return NULL;
  }

  rewind(f);
  while (fscanf(f, "%lld %lld", &u_ll, &v_ll) == 2) {
    int u = (int)u_ll, v = (int)v_ll;
    if ((flags & GRAPH_FLAG_IGNORE_SELFLOOPS) && u == v)
      continue;
    degrees[u]++;
    if (flags & GRAPH_FLAG_UNDIRECTED && u != v)
      degrees[v]++;
  }

  g->offsets[0] = 0;
  for (int i = 0; i < num_nodes; i++) {
    if ((long long)g->offsets[i] > INT_MAX - (long long)degrees[i]) {
      set_error("Offset overflow");
      free(degrees);
      graph_free(g);
      fclose(f);
      return NULL;
    }
    g->offsets[i + 1] = g->offsets[i] + (int)degrees[i];
  }

  g->m = g->offsets[num_nodes];
  g->edges = (int *)malloc(g->m * sizeof(int));
  if (!g->edges && g->m > 0) {
    set_error("Failed to allocate edges");
    free(degrees);
    graph_free(g);
    fclose(f);
    return NULL;
  }

  int *cursor = (int *)malloc(num_nodes * sizeof(int));
  if (!cursor) {
    set_error("Allocation failure (cursor)");
    free(degrees);
    graph_free(g);
    fclose(f);
    return NULL;
  }
  memcpy(cursor, g->offsets, num_nodes * sizeof(int));

  rewind(f);
  while (fscanf(f, "%lld %lld", &u_ll, &v_ll) == 2) {
    int u = (int)u_ll, v = (int)v_ll;
    if ((flags & GRAPH_FLAG_IGNORE_SELFLOOPS) && u == v)
      continue;
    g->edges[cursor[u]++] = v;
    if (flags & GRAPH_FLAG_UNDIRECTED && u != v) {
      g->edges[cursor[v]++] = u;
    }
  }

  free(degrees);
  free(cursor);
  fclose(f);

  if (flags & GRAPH_FLAG_DEDUP_EDGES) {
    int global_write_pos = 0;
    int read_start = g->offsets[0];

    for (int u = 0; u < num_nodes; u++) {
      int read_end = g->offsets[u + 1];
      int len = read_end - read_start;
      int start_idx = read_start;

      g->offsets[u] = global_write_pos;

      if (len > 0) {
        qsort(&g->edges[start_idx], len, sizeof(int), int_cmp);

        int w = 0;
        if (len > 0) {
          g->edges[start_idx + w++] = g->edges[start_idx];
          for (int r = 1; r < len; r++) {
            if (g->edges[start_idx + r] != g->edges[start_idx + w - 1]) {
              g->edges[start_idx + w++] = g->edges[start_idx + r];
            }
          }
        }
        int unique_len = w;

        if (global_write_pos != start_idx) {
          memmove(&g->edges[global_write_pos], &g->edges[start_idx],
                  unique_len * sizeof(int));
        }

        global_write_pos += unique_len;
      }

      read_start = read_end;
    }

    g->offsets[num_nodes] = global_write_pos;
    g->m = global_write_pos;
  }

  return g;
}

void graph_set_metadata_free_fn(Graph *g, graph_metadata_free_fn fn) {
  if (g)
    g->meta_free = fn;
}

void graph_free(Graph *g) {
  if (!g)
    return;

  if (g->nodes) {
    if (g->meta_free) {
      for (int i = 0; i < g->n; i++) {
        if (g->nodes[i].metadata)
          g->meta_free(g->nodes[i].metadata);
      }
    }
    free(g->nodes);
  }
  if (g->offsets)
    free(g->offsets);
  if (g->edges)
    free(g->edges);
  free(g);
}

int graph_num_nodes(const Graph *g) { return g ? g->n : 0; }

int graph_num_edges(const Graph *g) { return g ? g->m : 0; }

int graph_degree(const Graph *g, int node) {
  if (!g || node < 0 || node >= g->n)
    return 0;
  return g->offsets[node + 1] - g->offsets[node];
}

const int *graph_neighbors(const Graph *g, int node, int *out_count) {
  if (!g || node < 0 || node >= g->n) {
    if (out_count)
      *out_count = 0;
    return NULL;
  }
  int start = g->offsets[node];
  int end = g->offsets[node + 1];
  if (out_count)
    *out_count = end - start;
  return &g->edges[start];
}

void *graph_node_metadata(const Graph *g, int node) {
  if (!g || node < 0 || node >= g->n)
    return NULL;
  return g->nodes[node].metadata;
}

void graph_set_node_metadata(Graph *g, int node, void *data) {
  if (!g || node < 0 || node >= g->n)
    return;
  g->nodes[node].metadata = data;
}

int graph_has_edge(const Graph *g, int u, int v) {
  if (!g || u < 0 || u >= g->n)
    return 0;
  int count;
  const int *nbrs = graph_neighbors(g, u, &count);
  for (int i = 0; i < count; i++) {
    if (nbrs[i] == v)
      return 1;
  }
  return 0;
}
