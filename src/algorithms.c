#include "graph.h"
#include <stdlib.h>

typedef struct {
  int parent;
  int rank;
} UnionFind;

static void union_find_init(UnionFind *uf, int n) {
  for (int i = 0; i < n; i++) {
    uf[i].parent = i;
    uf[i].rank = 0;
  }
}

static int union_find_find(UnionFind *uf, int x) {
  if (uf[x].parent != x) {
    uf[x].parent = union_find_find(uf, uf[x].parent);
  }
  return uf[x].parent;
}

static void union_find_union(UnionFind *uf, int x, int y) {
  int px = union_find_find(uf, x);
  int py = union_find_find(uf, y);
  if (px != py) {
    if (uf[px].rank < uf[py].rank) {
      uf[px].parent = py;
    } else if (uf[px].rank > uf[py].rank) {
      uf[py].parent = px;
    } else {
      uf[py].parent = px;
      uf[px].rank++;
    }
  }
}

int graph_connected_components(const Graph *g, int *component_ids) {
  if (!g || !component_ids)
    return 0;

  UnionFind *uf = (UnionFind *)malloc(g->n * sizeof(UnionFind));
  if (!uf)
    return 0;

  union_find_init(uf, g->n);

  for (int u = 0; u < g->n; u++) {
    int count;
    const int *nbrs = graph_neighbors(g, u, &count);
    for (int i = 0; i < count; i++) {
      union_find_union(uf, u, nbrs[i]);
    }
  }

  int num_components = 0;
  int *comp_map = (int *)malloc(g->n * sizeof(int));
  if (!comp_map) {
    free(uf);
    return 0;
  }

  for (int i = 0; i < g->n; i++) {
    comp_map[i] = -1;
  }

  for (int i = 0; i < g->n; i++) {
    int root = union_find_find(uf, i);
    if (comp_map[root] == -1) {
      comp_map[root] = num_components++;
    }
    component_ids[i] = comp_map[root];
  }

  free(uf);
  free(comp_map);

  return num_components;
}

void graph_degree_stats(const Graph *g, double *avg, int *min, int *max) {
  if (!g)
    return;

  int min_deg = g->n > 0 ? graph_degree(g, 0) : 0;
  int max_deg = min_deg;
  long long sum = 0;

  for (int i = 0; i < g->n; i++) {
    int deg = graph_degree(g, i);
    sum += deg;
    if (deg < min_deg)
      min_deg = deg;
    if (deg > max_deg)
      max_deg = deg;
  }

  if (avg)
    *avg = g->n > 0 ? (double)sum / g->n : 0.0;
  if (min)
    *min = min_deg;
  if (max)
    *max = max_deg;
}
