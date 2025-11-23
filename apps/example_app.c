// Example application demonstrating most of Graphene's API
#include "graph.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

// Simple BFS visitor
static void print_bfs_node(int node, int depth, void *ctx) {
  (void)ctx;
  printf("  visited node %d (depth %d)\n", node, depth);
}

// Simple DFS visitor
static void print_dfs_node(int node, void *ctx) {
  (void)ctx;
  printf("  dfs visit node %d\n", node);
}

// Example metadata structure
typedef struct {
  double weight;
  int mark;
} Meta;

static void metadata_free(void *data) { free(data); }

static void attach_metadata(Graph *g) {
  for (int i = 0; i < graph_num_nodes(g); i++) {
    Meta *m = malloc(sizeof(Meta));
    m->weight = (double)i * 0.5;
    m->mark = i % 2;
    graph_set_node_metadata(g, i, m);
  }
}

static void print_metadata(Graph *g) {
  for (int i = 0; i < graph_num_nodes(g); i++) {
    Meta *m = (Meta *)graph_node_metadata(g, i);
    if (m) {
      printf("  node %d metadata: weight=%.2f mark=%d\n", i, m->weight,
             m->mark);
    }
  }
}

static void print_neighbors(Graph *g, int node) {
  int count = 0;
  const int *nbrs = graph_neighbors(g, node, &count);
  printf("  neighbors of %d (%d): ", node, count);
  for (int i = 0; i < count; i++)
    printf("%d ", nbrs[i]);
  printf("\n");
}

static void show_degree_stats(Graph *g) {
  double avg;
  int min, max;
  graph_degree_stats(g, &avg, &min, &max);
  printf("  degree stats: avg=%.2f min=%d max=%d\n", avg, min, max);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Usage: %s <graph.edgelist>\n", argv[0]);
    return 1;
  }

  const char *path = argv[1];

  printf("\n=== Load directed graph ===\n");
  Graph *g = graph_load(path);
  if (!g) {
    fprintf(stderr, "Failed to load: %s\n", graph_last_error());
    return 1;
  }

  printf("Graph loaded. Nodes=%d Edges=%d\n", graph_num_nodes(g),
         graph_num_edges(g));

  print_neighbors(g, 0);
  show_degree_stats(g);

  // BFS and DFS
  printf("\nBFS from node 0 (depth limit 2):\n");
  graph_bfs(g, 0, 2, print_bfs_node, NULL);

  printf("\nDFS from node 0:\n");
  graph_dfs(g, 0, print_dfs_node, NULL);

  // Connected components
  int n = graph_num_nodes(g);
  int *comp = malloc(n * sizeof(int));
  int num_comp = graph_connected_components(g, comp);
  printf("\nConnected components: %d\n", num_comp);
  for (int i = 0; i < n; i++) {
    printf("  node %d -> component %d\n", i, comp[i]);
  }
  free(comp);

  // Metadata test
  printf("\nAttach metadata to nodes:\n");
  graph_set_metadata_free_fn(g, metadata_free);
  attach_metadata(g);
  print_metadata(g);

  graph_free(g);

  printf("\n=== Load undirected graph (ignore self loops + dedup) ===\n");
  unsigned int flags = GRAPH_FLAG_UNDIRECTED | GRAPH_FLAG_IGNORE_SELFLOOPS |
                       GRAPH_FLAG_DEDUP_EDGES;
  Graph *g2 = graph_load_with_flags(path, flags);
  if (!g2) {
    fprintf(stderr, "Failed to load undirected: %s\n", graph_last_error());
    return 1;
  }
  printf("Graph loaded (undirected). Nodes=%d Edges=%d\n", graph_num_nodes(g2),
         graph_num_edges(g2));

  print_neighbors(g2, 0);
  show_degree_stats(g2);

  // RNG utilities demo
  printf("\nRandom node sampling (seed=123): ");
  utils_seed_rng(123);
  for (int i = 0; i < 5; i++) {
    int node = utils_rand_int(0, graph_num_nodes(g2) - 1);
    printf("%d ", node);
  }
  printf("\n");

  graph_free(g2);
  printf("\nDone.\n");
  return 0;
}
