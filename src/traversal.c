#include "graph.h"
#include <stdlib.h>
#include <string.h>

void graph_bfs(const Graph* g, int start, int max_depth,
               graph_visit_fn visit, void* ctx) {
    if (!g || start < 0 || start >= g->n || !visit) return;
    
    int* visited = (int*)calloc(g->n, sizeof(int));
    int* depths = (int*)calloc(g->n, sizeof(int));
    int* queue = (int*)malloc(g->n * sizeof(int));

    if (!visited || !depths || !queue) {
        if (visited) free(visited);
        if (depths) free(depths);
        if (queue) free(queue);
        return;
    }

    int head = 0, tail = 0;
    
    queue[tail++] = start;
    visited[start] = 1;
    depths[start] = 0;
    
    while (head < tail) {
        int u = queue[head++];
        int depth = depths[u];
        
        visit(u, depth, ctx);
        
        if (max_depth >= 0 && depth >= max_depth) continue;
        
        int count;
        const int* nbrs = graph_neighbors(g, u, &count);
        
        for (int i = 0; i < count; i++) {
            int v = nbrs[i];
            if (!visited[v]) {
                visited[v] = 1;
                depths[v] = depth + 1;
                queue[tail++] = v;
            }
        }
    }
    
    free(visited);
    free(depths);
    free(queue);
}

static void dfs_recursive(const Graph* g, int node, int* visited,
                         graph_visit_simple_fn visit, void* ctx) {
    visited[node] = 1;
    visit(node, ctx);
    
    int count;
    const int* nbrs = graph_neighbors(g, node, &count);
    
    for (int i = 0; i < count; i++) {
        int v = nbrs[i];
        if (!visited[v]) {
            dfs_recursive(g, v, visited, visit, ctx);
        }
    }
}

void graph_dfs(const Graph* g, int start,
               graph_visit_simple_fn visit, void* ctx) {
    if (!g || start < 0 || start >= g->n || !visit) return;
    
    int* visited = (int*)calloc(g->n, sizeof(int));
    if (!visited) return;

    dfs_recursive(g, start, visited, visit, ctx);
    free(visited);
}

