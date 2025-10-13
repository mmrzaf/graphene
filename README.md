# Graphene — Lightweight Graph Graph Library

**Graphene** is a lightweight, robust C library for graph representation and manipulation. Its design prioritizes **simplicity, modularity, and reproducibility**, while providing slightly richer functionality to support research, experiments, and algorithm development.

The core library provides **graph storage, traversal, node metadata, and small utilities**. Applications and research algorithms, such as metaheuristic influence maximization, sit separately in `apps/` and interact with the core through a stable API.

---

## Features

* **CSR Graph Core:** Memory-efficient representation for fast traversal.
* **Node Metadata:** Attach arbitrary payloads to nodes for experiments.
* **Graph Primitives:** Degree, neighbors, edge existence, connected components.
* **Traversal:** BFS with depth, DFS.
* **Utilities & Small Algorithms:** Degree statistics, frontier computations, connected components.
* **Modular Design:** Core library is decoupled from applications.
* **Reproducibility:** Deterministic RNG and seedset utilities for experiments.

---

**Notes:**

* `src/` is the core library; unaware of any application logic.
* `apps/` contains experiments or research algorithms, each with its own `main.c`.
* `examples/` provides sample graphs and seedsets.
* `tests/` validates correctness of graph operations.

---

## Core Graph API

### Graph Structure

```c
typedef struct Node {
    int id;
    void *metadata;   // Flexible payload per node
} Node;

typedef struct Graph {
    int n;            // Number of nodes
    int m;            // Number of edges
    int *offsets;     // CSR row offsets
    int *edges;       // Flattened adjacency list
    Node *nodes;      // Node metadata array
} Graph;
```

---

### Construction / Cleanup

```c
Graph* graph_load(const char* path);   // Load from edgelist
void graph_free(Graph* g);             // Free memory
```

---

### Node Queries

```c
int graph_num_nodes(const Graph* g);
int graph_num_edges(const Graph* g);
int graph_degree(const Graph* g, int node);
const int* graph_neighbors(const Graph* g, int node, int* out_count);
void* graph_node_metadata(const Graph* g, int node);
void graph_set_node_metadata(Graph* g, int node, void* data);
```

---

### Traversal

```c
// BFS with depth
void graph_bfs(const Graph* g, int start, int max_depth,
               void (*visit)(int node, int depth, void* ctx),
               void* ctx);

// DFS
void graph_dfs(const Graph* g, int start,
               void (*visit)(int node, void* ctx),
               void* ctx);
```

---

### Utilities / Small Algorithms

```c
int graph_has_edge(const Graph* g, int u, int v);
int graph_connected_components(const Graph* g, int* component_ids);
void graph_degree_stats(const Graph* g, double* avg, int* min, int* max);
```

**Optional helpers:**

* BFS-based distance counts
* Frontier computation for influence experiments
* Degree or centrality approximations

These small utilities provide ready-to-use functionality for research without bloating the core.

---

## Building Graphene

Simple Makefile-based build:

```bash
git clone https://github.com/darius-tg/graphene
cd graphene
make           # Builds libgraphene.a
```

Applications link against the core library:

```bash
cd apps/gwo_influence
make
./gwo_influence --graph ../../examples/small.edgelist --k 5 --pop 20 --iters 50
```

---

## Adding Applications

* Place research experiments or algorithms in `apps/`.
* Each app has its own `main.c` and may maintain its own Makefile.
* Apps use the core API (`graph.h` and `utils.h`) to access graph data and traversal functions.
* This separation keeps the library small and auditable.

---

## Examples

* `examples/small.edgelist` — tiny sample graph (0-based IDs)
* `examples/init_seedset.txt` — initial seedsets for deterministic experiments

**Sample usage:**

```c
#include "graph.h"

Graph *g = graph_load("examples/small.edgelist");

int n = graph_num_nodes(g);
printf("Graph has %d nodes\n", n);

int count;
const int *nbrs = graph_neighbors(g, 0, &count);
printf("Node 0 has %d neighbors\n", count);

graph_bfs(g, 0, 2, [](int node, int depth, void* ctx){
    printf("Visited %d at depth %d\n", node, depth);
}, NULL);

graph_free(g);
```

---

## Philosophy

* **Minimal & robust:** Core library does one thing — graph representation and traversal.
* **Slightly richer:** Node metadata, traversal with depth, connected components, small algorithms.
* **Separation of concerns:** Experiments live in `apps/`; core remains simple and reusable.
* **Extensible & research-ready:** Core provides all basic graph operations; research applications can build on top.

---

## Testing

* Located in `tests/`
* Validates CSR loading, BFS/DFS correctness, neighbor iteration, degree counts.
* Run tests:

```bash
make test
```

---

## License

MIT License — free for research, teaching, and personal use.

* `examples/` small graph and seedset
* `Makefile` and minimal tests

Do you want me to produce that starter repo now?
