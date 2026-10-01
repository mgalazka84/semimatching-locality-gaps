"""Independent checks on actual graphs, without using potential or state transitions.

Enumerate every labelled one-parent digraph on at most six servers by default.
For each server v, client v may use v or parent[v]; its current assignment is
parent[v] and assigning client v to v certifies OPT = number of clients.
Directed reachability determines the shortest improving alternating path.
This is a finite regression check, not a proof for arbitrary graph size.
"""
from collections import deque
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify_general_gap import gap, layered, materialize  # noqa: E402


@lru_cache(maxsize=None)
def exact_bound(r, K):
    return gap(r, K)


def actual_graph(parent):
    n = len(parent)
    children = [[] for _ in parent]
    for client, source in enumerate(parent):
        children[source].append(client)
    loads = [len(neighbours) for neighbours in children]
    cost = sum(a * (a + 1) // 2 for a in loads)
    shortest = None
    for source, a in enumerate(loads):
        if a < 2:
            continue
        distances = {source: 0}
        queue = deque([source])
        while queue:
            v = queue.popleft()
            if loads[v] <= a - 2:
                distance = distances[v]
                shortest = distance if shortest is None else min(shortest, distance)
                break
            for target in children[v]:
                if target not in distances:
                    distances[target] = distances[v] + 1
                    queue.append(target)
    return Fraction(cost, n), max(loads), shortest


def check_tree(eligibility, servers, require_path=False):
    """Check connectedness and edge count directly in the bipartite graph."""
    clients = len(eligibility)
    adjacency = [[] for _ in range(servers + clients)]
    edges = 0
    for u, menu in enumerate(eligibility):
        assert len(menu) == len(set(menu))
        for v in menu:
            assert 0 <= v < servers
            adjacency[servers + u].append(v)
            adjacency[v].append(servers + u)
            edges += 1
    seen = {0}
    queue = deque([0])
    while queue:
        for v in adjacency[queue.popleft()]:
            if v not in seen:
                seen.add(v)
                queue.append(v)
    assert len(seen) == servers + clients
    assert edges == servers + clients - 1
    if require_path:
        assert max(map(len, adjacency)) <= 2
    return edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-servers', type=int, default=6,
                        help='enumeration limit, 1 through 7 (default: 6)')
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results')
    args = parser.parse_args()
    if not __debug__:
        parser.error('verification requires assertions; do not use python -O or PYTHONOPTIMIZE')
    if not 1 <= args.max_servers <= 7:
        parser.error('--max-servers must be between 1 and 7')
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    enumerated = 0
    stable_checks = 0
    maxima = {}
    counts_by_size = []
    for n in range(1, args.max_servers + 1):
        count = 0
        for parent in product(range(n), repeat=n):
            ratio, K, shortest = actual_graph(parent)
            enumerated += 1
            count += 1
            for r in range(1, 5):
                if shortest is not None and shortest <= r:
                    continue
                stable_checks += 1
                assert ratio <= exact_bound(r, K), (parent, r, K, ratio)
                key = (n, r, K)
                maxima[key] = max(maxima.get(key, Fraction(0)), ratio)
        counts_by_size.append(dict(servers=n, graphs=count))

    # A different sharpness construction: an actual bipartite path with r+2
    # servers, one fixed client at the first server, and one client per link.
    path_checks = []
    for r in range(1, 13):
        n = r + 2
        parent = [0] + list(range(n - 1))
        eligibility = [sorted({v, parent[v]}) for v in range(n)]
        ratio, K, shortest = actual_graph(parent)
        assert K == 2 and shortest == r + 1
        assert ratio == Fraction(r + 3, r + 2) == exact_bound(r, 2)
        edges = check_tree(eligibility, n, require_path=True)
        path_checks.append(dict(r=r, servers=n, clients=n, edges=edges,
                                ratio=str(ratio), shortest_improving_path=shortest))
        if r == 3:
            witness = dict(r=r, K=2, servers=n, clients=n, eligibility=eligibility,
                           current_assignment=parent, optimal_assignment=list(range(n)),
                           cost=r+3, optimum=n, ratio=str(ratio))
            (output / 'r3_path_witness_5_servers.json').write_text(
                json.dumps(witness, indent=2) + '\n', encoding='utf-8')

    # Inspect the real bipartite edge sets of the materialized sharp examples.
    tree_checks = []
    for r, K in [(1, 2), (2, 3), (3, 3), (4, 2), (7, 2)]:
        graph = materialize(layered(r, K))
        edges = check_tree(graph['eligibility'], graph['servers'])
        server_degrees = [0] * graph['servers']
        for menu in graph['eligibility']:
            assert len(menu) <= 2
            for server in menu:
                server_degrees[server] += 1
        assert max(server_degrees) <= K
        tree_checks.append(dict(r=r, K=K, servers=graph['servers'],
                                clients=graph['clients'], edges=edges, is_tree=True,
                                maximum_server_degree=max(server_degrees)))

    report = dict(status='PASS', max_servers=args.max_servers,
                  enumerated_graphs=enumerated, stable_graph_parameter_checks=stable_checks,
                  counts_by_size=counts_by_size,
                  exhaustive_maxima=[dict(servers=n, r=r, K=K, observed_ratio=str(ratio),
                                          upper_bound=str(exact_bound(r, K)))
                                     for (n, r, K), ratio in sorted(maxima.items())],
                  path_witness_checks=path_checks, layered_tree_checks=tree_checks)
    (output / 'small_graph_verification.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f'PASS: {enumerated} actual parent-map graphs; {stable_checks} stable '
          f'graph/parameter checks; {len(path_checks)} sharp path witnesses; '
          f'{len(tree_checks)} layered bipartite trees.')
    print(f'Wrote independent graph report and path witness to {output}')


if __name__ == '__main__':
    main()
