"""Exact verification for matchable r-non-swappable semi-matchings.

Default checks use only the Python standard library. With --lp, SciPy checks
an independently assembled finite-state relaxation. Computations supplement
the all-r proof; they do not establish publication priority.
"""
from collections import deque
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from math import factorial, lcm
from pathlib import Path
import argparse
import csv
import json

DEFAULT_OUTPUT = Path(__file__).resolve().parent / 'results'


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer')


def _parameters(r, K, minimum_K=1):
    _positive_integer(r, 'r')
    _positive_integer(K, 'K')
    if K < minimum_K:
        raise ValueError(f'K must be at least {minimum_K}')


def gap(r, K):
    _parameters(r, K)
    if K == 1:
        return Q(1)
    E = sum((Q(1, factorial(a) ** r) for a in range(2, K)), Q(0))
    T = sum((Q(1, a * (a - 1) * factorial(a) ** r) for a in range(2, K)), Q(0))
    return 1 + (2 + E) / (2 * (r + 2 - T))


def states(r, K):
    _parameters(r, K)
    return [(0, r - 1)] + [(a, j) for a in range(1, K + 1) for j in range(r)]


def successors(r, K, state):
    _parameters(r, K)
    if state not in states(r, K):
        raise ValueError('state must belong to states(r, K)')
    a, j = state
    if a == 0:
        return []
    if j == 0:
        return [(a - 1, r - 1)] + [(b, 0) for b in range(a, K + 1)]
    return [(a, j - 1)] + [(b, 0) for b in range(a + 1, K + 1)]


def potential(r, K):
    _parameters(r, K, minimum_K=2)
    R = gap(r, K)
    eps = R - 1
    P = {(0, r - 1): -R}
    for j in range(r):
        P[1, j] = -1 - (j + 2) * eps
    for a in range(2, K + 1):
        s = -Q(a, 2) - 1 + eps / (a - 1)
        scale = factorial(a) * factorial(a - 1) ** (r - 1)
        delta = scale * sum(((Q(1, 2) + eps / (k * (k - 1))) /
                             factorial(k) ** r for k in range(a, K)), Q(0))
        assert delta >= 0
        for j in range(r):
            P[a, j] = s - a ** j * delta
    return R, P


def check_certificate(r, K):
    R, P = potential(r, K)
    slacks = []
    for state in states(r, K):
        a, j = state
        cost = Q(a * (a + 1), 2)
        if a == 0:
            assert P[state] == -R
        for target in successors(r, K, state):
            slack = P[state] - cost + R - a * P[target]
            assert slack >= 0, (r, K, state, target, slack)
            slacks.append(slack)
    return dict(r=r, K=K, exact_gap=str(R), states=len(P),
                tested_transitions=len(slacks), minimum_slack=str(min(slacks)), status='PASS')


def layered(r, K):
    _parameters(r, K, minimum_K=2)
    # Keep r layers at loads 1 through K-1 and just one top layer at K.
    L = r * (K - 1) + 1
    loads = [0] + [(i + r - 1) // r for i in range(1, L + 1)]
    raw = [Q(1)]
    for i in range(1, L):
        raw.append(raw[-1] / loads[i])
    raw.append(raw[-1] / (K - 1))
    scale = lcm(*(x.denominator for x in raw))
    counts = [int(scale * x) for x in raw]
    assert counts[-1] == 1
    n = sum(a * count for a, count in zip(loads, counts))
    m = sum(counts)
    cost = sum(a * (a + 1) // 2 * count for a, count in zip(loads, counts))
    assert n == m
    assert Q(cost, n) == gap(r, K)
    for i in range(1, L):
        assert loads[i] * counts[i] == counts[i - 1]
    assert (K - 1) * counts[L] == counts[L - 1]
    for i in range(L + 1):
        for t in range(1, min(r, i) + 1):
            assert loads[i] - loads[i - t] <= 1
    return dict(r=r, K=K, layers=L + 1, loads=loads, counts=counts,
                servers=m, clients=n, cost=cost, optimum=n, ratio=str(Q(cost, n)))


def materialize(data):
    levels, offset = [], 0
    for count in data['counts']:
        levels.append(list(range(offset, offset + count)))
        offset += count
    menus, current, comparison = [], [], []
    for i in range(1, len(levels)):
        number = data['loads'][i] - (i == len(levels) - 1)
        source = [v for v in levels[i] for _ in range(number)]
        assert len(source) == len(levels[i - 1])
        for v, w in zip(source, levels[i - 1]):
            menus.append([v, w]); current.append(v); comparison.append(w)
    for v in levels[-1]:
        menus.append([v]); current.append(v); comparison.append(v)
    return dict(**data, eligibility=menus, current_assignment=current,
                optimal_assignment=comparison)


def check_graph(graph):
    n, m, r, K = (graph[x] for x in ['clients', 'servers', 'r', 'K'])
    current, comparison = graph['current_assignment'], graph['optimal_assignment']
    loads, optloads, degrees = [0] * m, [0] * m, [0] * m
    parent, arcs = [None] * m, [[] for _ in range(m)]
    assert len(current) == len(comparison) == len(graph['eligibility']) == n
    for u, menu in enumerate(graph['eligibility']):
        v, w = current[u], comparison[u]
        assert v in menu and w in menu
        loads[v] += 1; optloads[w] += 1
        assert parent[w] is None
        parent[w] = v
        arcs[v].append(w)
        for z in menu:
            degrees[z] += 1
    assert optloads == [1] * m
    assert sum(a * (a + 1) // 2 for a in loads) == graph['cost']
    shortest = None
    for start in range(m):
        dist, queue = {start: 0}, deque([start])
        while queue:
            v = queue.popleft()
            if loads[start] >= loads[v] + 2:
                shortest = dist[v] if shortest is None else min(shortest, dist[v])
            for w in arcs[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; queue.append(w)
    assert shortest == r + 1
    # Recover states from actual ancestors, including the top self-loop.
    labels = []
    for v in range(m):
        ancestor, nearest = v, None
        for distance in range(1, r):
            ancestor = parent[ancestor]
            assert loads[ancestor] <= loads[v] + 1
            if nearest is None and loads[ancestor] == loads[v] + 1:
                nearest = distance
        labels.append((loads[v], 0 if nearest is None else r - nearest))
    R, P = potential(r, K)
    for v in range(m):
        allowed = successors(r, K, labels[v])
        for w in arcs[v]:
            assert labels[w] in allowed
        residual = P[labels[v]] - Q(loads[v] * (loads[v] + 1), 2) + R
        residual -= sum((P[labels[w]] for w in arcs[v]), Q(0))
        assert residual == 0  # Every vertex of this extremal witness is tight.
    return dict(r=r, K=K, clients=n, servers=m, cost=graph['cost'], optimum=n,
                ratio=graph['ratio'], shortest_improving_path=shortest,
                maximum_server_degree=max(degrees),
                maximum_client_degree=max(map(len, graph['eligibility'])),
                recovered_states=sorted(set(labels)), all_vertex_certificate_slacks='0', status='PASS')


def infinite_interval(r, digits=30):
    """Rigorous rational enclosure, using a geometric bound for the series tail."""
    _positive_integer(r, 'r')
    _positive_integer(digits, 'digits')
    tolerance = Q(1, 10 ** digits)
    L = 2
    while True:
        E = sum((Q(1, factorial(a) ** r) for a in range(2, L + 1)), Q(0))
        T = sum((Q(1, a * (a - 1) * factorial(a) ** r)
                 for a in range(2, L + 1)), Q(0))
        tail_E = Q(1, factorial(L + 1) ** r) / (1 - Q(1, (L + 2) ** r))
        tail_T = tail_E / (L * (L + 1))
        lower = 1 + (2 + E) / (2 * (r + 2 - T))
        upper = 1 + (2 + E + tail_E) / (2 * (r + 2 - T - tail_T))
        assert lower < upper
        if upper - lower < tolerance:
            break
        L += 1
    if r == 1:
        assert lower < Q(3, 2) < upper
    with localcontext() as ctx:
        ctx.prec = digits + 10
        midpoint = (lower + upper) / 2
        decimal = Decimal(midpoint.numerator) / Decimal(midpoint.denominator)
        display = f'{decimal:.15f}'
    return dict(r=r, decimal_15_places=display, included_through=L,
                exact_lower=str(lower), exact_upper=str(upper), exact_width=str(upper-lower))


def numerical_lp(r, K):
    _parameters(r, K)
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix
    types = states(r, K)
    index = {state: i for i, state in enumerate(types)}
    edges = [(index[state], index[target]) for state in types
             for target in successors(r, K, state)]
    S, L = len(types), len(types) + len(edges)
    A = lil_matrix((2 * S + 1, L))
    for i, (a, j) in enumerate(types):
        A[i, i] = -a; A[S + i, i] = -1; A[2 * S, i] = 1
    for e, (source, target) in enumerate(edges):
        A[source, S + e] = 1; A[S + target, S + e] = 1
    solution = linprog([-a * (a + 1) / 2 for a, j in types] + [0.] * len(edges),
                       A_eq=A.tocsr(), b_eq=[0.] * (2 * S) + [1.],
                       bounds=(0, None), method='highs')
    assert solution.success, solution.message
    assert abs(-solution.fun - float(gap(r, K))) < 1e-8
    return dict(r=r, K=K, numerical_optimum=-float(solution.fun),
                exact_formula=str(gap(r, K)), status='PASS; floating-point corroboration only')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lp', action='store_true',
                        help='also run 24 floating-point LP checks (requires SciPy)')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT,
                        help='directory for generated JSON/CSV files (default: results/ beside script)')
    args = parser.parse_args()
    if not __debug__:
        parser.error('verification requires assertions; do not use python -O or PYTHONOPTIMIZE')
    if args.lp:
        try:
            import scipy  # noqa: F401
        except ImportError:
            parser.error('--lp requires SciPy; install requirements-optional.txt')
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = [(r, K) for r in range(1, 13) for K in range(2, 11)]
    certificates = [check_certificate(r, K) for r, K in cases]
    constructions = [layered(r, K) for r, K in cases]
    assert all(gap(r, 1) == 1 for r in range(1, 13))
    graph_checks = []
    for r, K in [(1, 2), (2, 3), (3, 3), (4, 2), (7, 2)]:
        graph = materialize(layered(r, K))
        graph_checks.append(check_graph(graph))
        if (r, K) == (3, 3):
            (output / 'r3_witness_79_servers.json').write_text(json.dumps(graph, indent=2) + '\n', encoding='utf-8')
    intervals = [infinite_interval(r) for r in [1, 2, 3, 4, 5, 8, 10, 20, 50, 100]]
    results = dict(exact_certificates=certificates, exact_construction_count=len(constructions),
                   explicit_graph_checks=graph_checks, infinite_constant_intervals=intervals)
    if args.lp:
        results['numerical_lp'] = [numerical_lp(r, K) for r in [1, 2, 3, 4, 5, 8]
                                   for K in [2, 3, 4, 6]]
    (output / 'general_gap_verification.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
    with (output / 'exact_gap_table.csv').open('w', newline='', encoding='utf-8') as handle:
        fields = ['r', 'K', 'servers', 'clients', 'cost', 'optimum', 'ratio']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: data[k] for k in fields} for data in constructions)
    with (output / 'infinite_gap_table.csv').open('w', newline='', encoding='utf-8') as handle:
        fields = ['r', 'decimal_15_places', 'included_through', 'exact_width']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: data[k] for k in fields} for data in intervals)
    print(f'PASS: {len(certificates)} exact certificates; '
          f'{sum(x["tested_transitions"] for x in certificates)} transition inequalities; '
          f'{len(constructions)} exact constructions; {len(graph_checks)} explicit graphs; '
          f'{len(intervals)} rigorous constant enclosures; '
          f'{len(results.get("numerical_lp", []))} numerical LP checks.')
    print(f'Wrote verification report, witness, and two tables to {output}')


if __name__ == '__main__':
    main()
