#include "graph.h"
#include "utils.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_passed = 0;
static int test_failed = 0;

#define TEST(name)                                                             \
  printf("Running test: %s...", name);                                         \
  fflush(stdout);

#define PASS()                                                                 \
  printf(" PASS\n");                                                           \
  test_passed++;

#define FAIL(msg)                                                              \
  printf(" FAIL: %s\n", msg);                                                  \
  test_failed++;

void test_graph_creation() {
  TEST("graph_creation");

  Graph *g = graph_create(10);
  assert(g != NULL);
  assert(graph_num_nodes(g) == 10);
  assert(graph_num_edges(g) == 0);

  graph_free(g);
  PASS();
}

void test_graph_load() {
  TEST("graph_load");

  // Create a temporary test file
  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n1 2\n2 3\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");
  assert(g != NULL);
  assert(graph_num_nodes(g) == 4);
  assert(graph_num_edges(g) == 4);

  assert(graph_degree(g, 0) == 2);
  assert(graph_degree(g, 1) == 1);
  assert(graph_degree(g, 2) == 1);
  assert(graph_degree(g, 3) == 0);
  // assert(graph_degree(g, 0) == 2);
  // assert(graph_degree(g, 1) == 3);
  // assert(graph_degree(g, 2) == 4);
  // assert(graph_degree(g, 3) == 3);

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

void test_neighbors() {
  TEST("neighbors");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n0 3\n1 2\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  int count;
  const int *nbrs = graph_neighbors(g, 0, &count);
  assert(count == 3);

  // Check neighbors of node 0
  int found[3] = {0, 0, 0};
  for (int i = 0; i < count; i++) {
    if (nbrs[i] == 1)
      found[0] = 1;
    if (nbrs[i] == 2)
      found[1] = 1;
    if (nbrs[i] == 3)
      found[2] = 1;
  }
  assert(found[0] && found[1] && found[2]);

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

void test_has_edge() {
  TEST("has_edge");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n1 3\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  assert(graph_has_edge(g, 0, 1) == 1);
  assert(graph_has_edge(g, 0, 2) == 1);
  assert(graph_has_edge(g, 1, 3) == 1);
  assert(graph_has_edge(g, 0, 3) == 0);
  assert(graph_has_edge(g, 2, 1) == 0);

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

static int bfs_visit_count = 0;
static void bfs_visit(int node, int depth, void *ctx) { bfs_visit_count++; }

void test_bfs() {
  TEST("bfs");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n1 3\n2 3\n3 4\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  bfs_visit_count = 0;
  graph_bfs(g, 0, -1, bfs_visit, NULL);
  assert(bfs_visit_count == 5);

  bfs_visit_count = 0;
  graph_bfs(g, 0, 1, bfs_visit, NULL);
  assert(bfs_visit_count == 3); // Node 0, 1, 2

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

static int dfs_visit_count = 0;
static void dfs_visit(int node, void *ctx) { dfs_visit_count++; }

void test_dfs() {
  TEST("dfs");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n1 3\n2 4\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  dfs_visit_count = 0;
  graph_dfs(g, 0, dfs_visit, NULL);
  assert(dfs_visit_count == 5);

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

void test_connected_components() {
  TEST("connected_components");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n1 2\n3 4\n5 6\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  int *components = (int *)malloc(g->n * sizeof(int));
  int num_comp = graph_connected_components(g, components);

  assert(num_comp == 3);
  assert(components[0] == components[1]);
  assert(components[1] == components[2]);
  assert(components[3] == components[4]);
  assert(components[5] == components[6]);

  free(components);
  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

void test_degree_stats() {
  TEST("degree_stats");

  FILE *f = fopen("test_temp.edgelist", "w");
  fprintf(f, "0 1\n0 2\n0 3\n1 2\n");
  fclose(f);

  Graph *g = graph_load("test_temp.edgelist");

  double avg;
  int min, max;
  graph_degree_stats(g, &avg, &min, &max);

  assert(min == 0);
  assert(max == 3);
  assert(avg > 0.9 && avg < 1.1); // ~1.0

  graph_free(g);
  remove("test_temp.edgelist");
  PASS();
}

void test_metadata() {
  TEST("metadata");

  Graph *g = graph_create(5);

  int data1 = 42;
  int data2 = 99;

  graph_set_node_metadata(g, 0, &data1);
  graph_set_node_metadata(g, 1, &data2);

  int *retrieved1 = (int *)graph_node_metadata(g, 0);
  int *retrieved2 = (int *)graph_node_metadata(g, 1);

  assert(*retrieved1 == 42);
  assert(*retrieved2 == 99);
  assert(graph_node_metadata(g, 2) == NULL);

  graph_free(g);
  PASS();
}

void test_utils_rng() {
  TEST("utils_rng");

  utils_seed_rng(12345);
  int r1 = utils_rand_int(0, 100);

  utils_seed_rng(12345);
  int r2 = utils_rand_int(0, 100);

  assert(r1 == r2); // Deterministic

  double d = utils_rand_double();
  assert(d >= 0.0 && d <= 1.0);

  PASS();
}

int main() {
  printf("Running Graphene Tests\n");
  printf("======================\n\n");

  test_graph_creation();
  test_graph_load();
  test_neighbors();
  test_has_edge();
  test_bfs();
  test_dfs();
  test_connected_components();
  test_degree_stats();
  test_metadata();
  test_utils_rng();

  printf("\n======================\n");
  printf("Results: %d passed, %d failed\n", test_passed, test_failed);

  return test_failed > 0 ? 1 : 0;
}
