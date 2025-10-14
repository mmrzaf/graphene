#ifndef GRAPHENE_GRAPH_H
#define GRAPHENE_GRAPH_H

#include <stddef.h>

// Graph load flags
#define GRAPH_FLAG_NONE              0u
#define GRAPH_FLAG_UNDIRECTED        1u
#define GRAPH_FLAG_IGNORE_SELFLOOPS  2u
#define GRAPH_FLAG_DEDUP_EDGES       4u

typedef struct Node {
    int id;
    void *metadata;
} Node;

typedef void (*graph_metadata_free_fn)(void*);

typedef struct Graph {
    int n;                     // number of nodes
    int m;                     // number of edges
    int *offsets;              // CSR offsets
    int *edges;                // adjacency list
    Node *nodes;               // metadata per node
    unsigned int flags;        // graph flags
    graph_metadata_free_fn meta_free; // optional metadata destructor
} Graph;

// Construction / Cleanup
Graph* graph_create(int num_nodes);
Graph* graph_load(const char* path);
Graph* graph_load_with_flags(const char* path, unsigned int flags);
void graph_free(Graph* g);
void graph_set_metadata_free_fn(Graph* g, graph_metadata_free_fn fn);

// Node Queries
int graph_num_nodes(const Graph* g);
int graph_num_edges(const Graph* g);
int graph_degree(const Graph* g, int node);
const int* graph_neighbors(const Graph* g, int node, int* out_count);
void* graph_node_metadata(const Graph* g, int node);
void graph_set_node_metadata(Graph* g, int node, void* data);

// Edge Operations
int graph_has_edge(const Graph* g, int u, int v);

// Traversal
typedef void (*graph_visit_fn)(int node, int depth, void* ctx);
typedef void (*graph_visit_simple_fn)(int node, void* ctx);

void graph_bfs(const Graph* g, int start, int max_depth,
               graph_visit_fn visit, void* ctx);
void graph_dfs(const Graph* g, int start,
               graph_visit_simple_fn visit, void* ctx);

// Algorithms
int graph_connected_components(const Graph* g, int* component_ids);
void graph_degree_stats(const Graph* g, double* avg, int* min, int* max);

// Errors
const char* graph_last_error(void);

#endif // GRAPHENE_GRAPH_H

