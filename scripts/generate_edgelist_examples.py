import os
from generate_edgelist import stream_edges

configs = [
    {"nodes": 100, "edges": 200, "directed": False, "weighted": False, "allow_self_loops": False, "multi_edges": False},
    {"nodes": 100, "edges": 300, "directed": False, "weighted": False, "weight_range": (1,10)},
    {"nodes": 500, "edges": 1000, "directed": False, "weighted": False, "allow_self_loops": True, "weight_range": (5,20)},
    {"nodes": 500, "edges": 1500, "directed": False, "weighted": False, "multi_edges": True},
    {"nodes": 1000, "edges": 5000, "directed": False, "weighted": False},
    {"nodes": 1000, "edges": 5000, "directed": False, "weighted": False, "weight_range": (10,50)},
    {"nodes": 5000, "edges": 20000, "directed": False, "weighted": False, "allow_self_loops": True, "multi_edges": True},
    {"nodes": 5000, "edges": 20000, "directed": False, "weighted": False, "weight_range": (1,100)},
    {"nodes": 10000, "edges": 50000, "directed": False, "weighted": False},
    {"nodes": 10000, "edges": 100000, "directed": False, "weighted": False, "multi_edges": True, "weight_range": (0.1,10)},
]

for cfg in configs:
    parts = [
        f"nodes{cfg['nodes']}",
        f"edges{cfg['edges']}",
        "dir" if cfg.get("directed", False) else "undir",
        "weighted" if cfg.get("weighted", False) else "unweighted",
    ]
    if cfg.get("allow_self_loops", False):
        parts.append("selfloop")
    if cfg.get("multi_edges", False):
        parts.append("multi")
    filename = "_".join(parts) + ".edgelist"

    print(f"Generating {filename} ...")
    stream_edges(
        num_nodes=cfg['nodes'],
        num_edges=cfg['edges'],
        directed=cfg.get("directed", False),
        weighted=cfg.get("weighted", False),
        weight_range=cfg.get("weight_range", (1,10)),
        allow_self_loops=cfg.get("allow_self_loops", False),
        allow_multi_edges=cfg.get("multi_edges", False),
        output_file=filename,
    )

