# Graphene — Lightweight Graph Library

Graphene is a lightweight, robust C library for graph representation and manipulation. Its design prioritizes simplicity, modularity, and reproducibility, while providing slightly richer functionality to support research, experiments, and algorithm development.

## Features

* **CSR Graph Core**: Memory-efficient compressed sparse row representation
* **Node Metadata**: Attach arbitrary payloads to nodes
* **Graph Primitives**: Degree queries, neighbor iteration, edge existence
* **Traversal**: BFS with depth tracking, DFS
* **Algorithms**: Connected components, degree statistics
* **Utilities**: Deterministic RNG, seedset I/O
* **Modular Design**: Core library decoupled from applications

## Building

### Using Make

```bash
make                # Build library
make test          # Run tests
make clean         # Clean build files
make install       # Install to /usr/local
```

### Using CMake

```bash
mkdir build && cd build
cmake ..
make
make test
sudo make install
```

## Quick Start

```c
#include "graph.h"
#include "utils.h"

int main() {
    // Load graph from edgelist
    Graph *g = graph_load("examples/small.edgelist");
    
    printf("Nodes: %d, Edges: %d\n", 
           graph_num_nodes(g), graph_num_edges(g));
    
    // Query neighbors
    int count;
    const int *neighbors = graph_neighbors(g, 0, &count);
    printf("Node 0 has %d neighbors\n", count);
    
    // BFS traversal
    graph_bfs(g, 0, 2, my_visit_fn, NULL);
    
    // Cleanup
    graph_free(g);
    return 0;
}
```

## API Reference

### Graph Construction

```c
Graph* graph_create(int num_nodes);
Graph* graph_load(const char* path);
void graph_free(Graph* g);
```

### Queries

```c
int graph_num_nodes(const Graph* g);
int graph_num_edges(const Graph* g);
int graph_degree(const Graph* g, int node);
const int* graph_neighbors(const Graph* g, int node, int* out_count);
int graph_has_edge(const Graph* g, int u, int v);
```

### Metadata

```c
void* graph_node_metadata(const Graph* g, int node);
void graph_set_node_metadata(Graph* g, int node, void* data);
```

### Traversal

```c
void graph_bfs(const Graph* g, int start, int max_depth,
               void (*visit)(int node, int depth, void* ctx),
               void* ctx);

void graph_dfs(const Graph* g, int start,
               void (*visit)(int node, void* ctx),
               void* ctx);
```

### Algorithms

```c
int graph_connected_components(const Graph* g, int* component_ids);
void graph_degree_stats(const Graph* g, double* avg, int* min, int* max);
```

### Utilities

```c
void utils_seed_rng(unsigned int seed);
int utils_rand_int(int min, int max);
double utils_rand_double(void);
void utils_shuffle_int_array(int* arr, int n);
int* utils_load_seedset(const char* path, int* out_count);
void utils_save_seedset(const char* path, const int* seeds, int count);
```

## File Format

Edge list format (space or tab separated):
```
0 1
0 2
1 3
...
```

Node IDs should be 0-indexed integers.

## Examples

### BFS with Custom Visitor

```c
void my_visitor(int node, int depth, void* ctx) {
    printf("Node %d at depth %d\n", node, depth);
}

graph_bfs(g, 0, -1, my_visitor, NULL);
```

### Finding Connected Components

```c
int* components = malloc(g->n * sizeof(int));
int num_comp = graph_connected_components(g, components);
printf("Graph has %d components\n", num_comp);
```

### Using Node Metadata

```c
typedef struct {
    double score;
    int active;
} NodeData;

NodeData* data = malloc(sizeof(NodeData));
data->score = 0.75;
data->active = 1;

graph_set_node_metadata(g, node_id, data);

// Later...
NodeData* retrieved = graph_node_metadata(g, node_id);
printf("Score: %f\n", retrieved->score);
```

### Deterministic Random Seeds

```c
utils_seed_rng(42);  // Reproducible experiments
int* seeds = malloc(k * sizeof(int));
for (int i = 0; i < k; i++) {
    seeds[i] = utils_rand_int(0, g->n - 1);
}
utils_shuffle_int_array(seeds, k);
```

## Adding Applications

Place your research code in `apps/`:

```
apps/
└── my_algorithm/
    ├── main.c
    └── Makefile
```

Link against libgraphene:

```makefile
CC = gcc
CFLAGS = -Wall -O2 -I../../include
LDFLAGS = -L../../lib -lgraphene

my_app: main.c
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@
```

## Testing

The test suite validates:
- Graph loading and CSR construction
- Neighbor queries and edge detection
- BFS/DFS traversal correctness
- Connected component detection
- Degree statistics
- Metadata storage
- Deterministic RNG

Run tests:
```bash
make test
```

## Design Philosophy

- **Minimal Core**: Graph representation and traversal only
- **No Dependencies**: Pure C99, standard library only
- **Memory Efficient**: CSR format for large sparse graphs
- **Extensible**: Metadata system for custom node data
- **Reproducible**: Deterministic RNG for experiments
- **Separation**: Applications separate from core library

## Performance Notes

- CSR format: O(1) neighbor access, O(degree) iteration
- BFS/DFS: O(V + E) time complexity
- Connected components: O(V + E) with union-find
- Memory: ~(V + 2E) integers for graph structure

## License

MIT License - Free for research, teaching, and personal use.

