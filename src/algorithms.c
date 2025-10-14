#include "graph.h"
#include <stdlib.h>

static void union_find_init(int* parent, int n) {
    for (int i = 0; i < n; i++) {
        parent[i] = i;
    }
}

static int union_find_find(int* parent, int x) {
    if (parent[x] != x) {
        parent[x] = union_find_find(parent, parent[x]);
    }
    return parent[x];
}

static void union_find_union(int* parent, int x, int y) {
    int px = union_find_find(parent, x);
    int py = union_find_find(parent, y);
    if (px != py) {
        parent[px] = py;
    }
}

int graph_connected_components(const Graph* g, int* component_ids) {
    if (!g || !component_ids) return 0;
    
    int* parent = (int*)malloc(g->n * sizeof(int));
    union_find_init(parent, g->n);
    
    // Union all edges
    for (int u = 0; u < g->n; u++) {
        int count;
        const int* nbrs = graph_neighbors(g, u, &count);
        for (int i = 0; i < count; i++) {
            union_find_union(parent, u, nbrs[i]);
        }
    }
    
    // Assign component IDs and count
    int num_components = 0;
    int* comp_map = (int*)malloc(g->n * sizeof(int));
    for (int i = 0; i < g->n; i++) {
        comp_map[i] = -1;
    }
    
    for (int i = 0; i < g->n; i++) {
        int root = union_find_find(parent, i);
        if (comp_map[root] == -1) {
            comp_map[root] = num_components++;
        }
        component_ids[i] = comp_map[root];
    }
    
    free(parent);
    free(comp_map);
    
    return num_components;
}

void graph_degree_stats(const Graph* g, double* avg, int* min, int* max) {
    if (!g) return;
    
    int min_deg = g->n > 0 ? graph_degree(g, 0) : 0;
    int max_deg = min_deg;
    long long sum = 0;
    
    for (int i = 0; i < g->n; i++) {
        int deg = graph_degree(g, i);
        sum += deg;
        if (deg < min_deg) min_deg = deg;
        if (deg > max_deg) max_deg = deg;
    }
    
    if (avg) *avg = g->n > 0 ? (double)sum / g->n : 0.0;
    if (min) *min = min_deg;
    if (max) *max = max_deg;
}

