import numpy as np
import argparse

def validate_params(num_nodes, num_edges, directed, allow_self_loops, allow_multi_edges):
    max_edges = num_nodes * (num_nodes - 1)
    if not directed:
        max_edges //= 2
    if allow_self_loops:
        max_edges += num_nodes
    if not allow_multi_edges and num_edges > max_edges:
        raise ValueError(
            f"Cannot generate {num_edges} unique edges with the current constraints. Max possible: {max_edges}"
        )

def stream_edges(
    num_nodes, num_edges,
    directed=False, weighted=False,
    weight_range=(1,10),
    allow_self_loops=False, allow_multi_edges=False,
    batch_size=10_000, seed=None,
    output_file="edgelist.txt"
):
    if seed is not None:
        np.random.seed(seed)

    validate_params(num_nodes, num_edges, directed, allow_self_loops, allow_multi_edges)

    edges_generated = 0
    unique_edges = set() if not allow_multi_edges else None

    with open(output_file, 'w') as f:
        while edges_generated < num_edges:
            current_batch = min(batch_size, num_edges - edges_generated)

            u = np.random.randint(0, num_nodes, size=current_batch)
            v = np.random.randint(0, num_nodes, size=current_batch)

            if not allow_self_loops:
                mask = u != v
                u, v = u[mask], v[mask]
                if len(u) == 0:
                    continue  # skip empty batch

            if not directed:
                uv = np.stack([u, v], axis=1)
                uv.sort(axis=1)
                u, v = uv[:,0], uv[:,1]

            for a, b in zip(u, v):
                edge = (a, b)
                if weighted:
                    w = np.random.uniform(weight_range[0], weight_range[1])
                    edge = (a, b, w)

                if not allow_multi_edges:
                    if edge in unique_edges:
                        continue
                    unique_edges.add(edge)

                f.write(' '.join(map(str, edge)) + '\n')
                edges_generated += 1
                if edges_generated >= num_edges:
                    break

            # Optional: clear set if memory grows too big for massive graphs
            if not allow_multi_edges and len(unique_edges) > 10_000_000:
                unique_edges.clear()

    print(f"Generated {edges_generated} edges -> {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Memory-efficient streaming edgelist generator")
    parser.add_argument("--nodes", type=int, default=1_000_000)
    parser.add_argument("--edges", type=int, default=10_000_000)
    parser.add_argument("--directed", action="store_true")
    parser.add_argument("--weighted", action="store_true")
    parser.add_argument("--weight_min", type=float, default=1)
    parser.add_argument("--weight_max", type=float, default=10)
    parser.add_argument("--allow_self_loops", action="store_true")
    parser.add_argument("--multi_edges", action="store_true")
    parser.add_argument("--batch_size", type=int, default=10_000)
    parser.add_argument("--output", type=str, default="edgelist.txt")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    stream_edges(
        num_nodes=args.nodes,
        num_edges=args.edges,
        directed=args.directed,
        weighted=args.weighted,
        weight_range=(args.weight_min, args.weight_max),
        allow_self_loops=args.allow_self_loops,
        allow_multi_edges=args.multi_edges,
        batch_size=args.batch_size,
        seed=args.seed,
        output_file=args.output
    )

