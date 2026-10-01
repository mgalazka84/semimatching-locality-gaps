# Exact locality gaps for matchable semi-matchings

Reproducible research code for a working exact characterization of the cost of
locally optimal semi-matchings. The code checks rational certificates, constructs
sharp examples, verifies actual alternating paths, and independently enumerates
small graphs. A finite computation supplements the mathematical proof; it does
not establish the theorem for all parameters or the publication priority of the
result.

## Model and result

A finite bipartite graph has clients and servers. Each client is assigned to one
eligible server. A server with current load `a` contributes `a(a+1)/2` to the
cost. The graph must admit a matching covering every client, so the optimum cost
is the number of clients `n`.

An assignment is **r-non-swappable** if no improving alternating path moves at
most `r` clients. Here `r` counts moved clients, so such a path has at most `2r`
bipartite edges. Moving along a path improves the cost exactly when the starting
server's load exceeds the ending server's load by at least two. `K` is an upper
bound on the **current assigned load**, not on graph degree.

For positive integers `r` and `K >= 2`, define

```math
E_{r,K}=\sum_{a=2}^{K-1}\frac{1}{(a!)^r},\qquad
T_{r,K}=\sum_{a=2}^{K-1}\frac{1}{a(a-1)(a!)^r}.
```

The exact worst ratio within this class is

```math
\Gamma^{\mathrm{match}}_{r,K}
=1+\frac{2+E_{r,K}}{2(r+2-T_{r,K})}.
```

Empty sums are zero. Thus `K=2` gives `1+1/(r+2)`; `K=1` gives `1`.
Every finite-`K` bound is attained by a finite bipartite **tree** with client
degree at most two and server degree at most `K`. For `K=2`, a bipartite path
with `r+2` clients and `r+2` servers suffices.

For unrestricted current loads, replacing both sums by their infinite versions
gives the supremum. The asymptotic excess over one is
`1/(r+2) + O(2^(-r)/r)` as `r` tends to infinity.

| r | Unrestricted-load supremum, rounded |
|---:|---:|
| 1 | 1.500000000000000 |
| 2 | 1.294503158976119 |
| 3 | 1.215700053924918 |
| 4 | 1.172843490377998 |
| 5 | 1.145423532052361 |
| 10 | 1.083377416781517 |

These statements concern client-matchable instances. The repository does not
claim the same formula for arbitrary semi-matching instances and does not
provide a new distributed algorithm. The exact characterization is a working
research result; its priority remains subject to literature review.

Related work includes Correa and Muñoz, *Performance guarantees of local
search for minsum scheduling problems*, Mathematical Programming **191**
(2022), 847–869, [doi:10.1007/s10107-020-01571-5](https://doi.org/10.1007/s10107-020-01571-5)
([author-hosted text](https://www.dii.uchile.cl/~jcorrea/papers/Journals/2022CM.pdf)).
That paper studies Jump and Swap local search. A Swap exchanging two jobs is a
different neighborhood from the alternating paths of up to `r` clients used
here; its results must be compared under their stated model assumptions.

## Run the verification

The exact checks require **Python 3.10 or later** and no external packages.
From the repository root:

```sh
python verify_general_gap.py
python tests/check_small_graphs.py
```

Both scripts write reports to the ignored `results/` directory next to the main
script. They leave the checked-in `data/` directory unchanged. Pass
`--output-dir PATH` to choose a different destination. Existing generated files
with the same names in that destination are overwritten. Run without `-O` or
`PYTHONOPTIMIZE`; both command-line entry points reject disabled assertions.

Optional floating-point linear-programming checks use SciPy. The pinned tested
version, SciPy 1.17.0, requires **Python 3.11 or later**:

```sh
python -m pip install -r requirements-optional.txt
python verify_general_gap.py --lp
```

To reproduce the checked-in reports, including optional LP results:

```sh
python verify_general_gap.py --lp --output-dir data
python tests/check_small_graphs.py --output-dir data
```

Floating-point LP objective values can differ in their last digits across
platforms. These checks use an absolute tolerance of `1e-8`; the main certificate
and construction checks use exact `fractions.Fraction` arithmetic.

## What is checked

`verify_general_gap.py` checks:

- 108 exact rational potential certificates: `r=1,...,12`, `K=2,...,10`.
- 17,730 inequalities for the allowed state transitions.
- 108 exact layered constructions attaining the finite formula.
- Five materialized bipartite graphs: valid assignments, a client-covering
  comparison matching, the exact cost, the shortest improving path, actual
  ancestor-derived states, and equality in every vertex certificate.
- Ten rigorous rational enclosures for the unrestricted-load constants, each
  narrower than `10^-30`, using geometric bounds on the omitted series tails.
- With `--lp`, 24 finite-state flow LP optima, compared with the exact formula.

`tests/check_small_graphs.py` checks the bound without using the potential or its
state-transition rules:

- All **50,069** labelled one-parent digraphs on one through six servers. Client
  `v` is assigned to `parent[v]` and may alternatively use server `v`, certifying
  that the optimum is `n`. Loops and directed cycles are included.
- Actual directed reachability determines the shortest improving alternating
  path. Every stable case for `r=1,...,4` is compared with the formula, giving
  **20,274** graph/parameter checks. This enumerates the stated graph family,
  not all bipartite eligibility graphs.
- Twelve independent sharp path examples, one for each `r=1,...,12`, including
  exact ratios, shortest improving paths, and the bipartite path property.
- Connectedness and tree edge counts of the five materialized layered examples.

The enumeration limit can be changed with `--max-servers N`, for `1 <= N <= 7`.
The number of parent maps at size `N` is `N^N`.

## Construction and files

The compressed layered construction uses `r` layers at each load from `1` to
`K-1`, one zero-load layer, and a single top server at load `K`. Each non-top
server has one comparison client coming from the next layer. The top server has
one fixed client; every other client has two eligible servers in adjacent
layers. The top server has `K-1` children. The resulting eligibility graph is a
tree.

| File | Contents |
|---|---|
| `verify_general_gap.py` | Exact formula, potentials, construction, witness checks, interval bounds, optional LP |
| `tests/check_small_graphs.py` | Exhaustive reachability checks and independent graph structure checks |
| `data/general_gap_verification.json` | Certificate, witness, interval, and LP verification results |
| `data/exact_gap_table.csv` | Exact finite gaps and construction sizes |
| `data/infinite_gap_table.csv` | Rounded infinite constants and rigorous interval widths |
| `data/r3_witness_79_servers.json` | Sharp `r=3, K=3` tree: 79 clients, 79 servers, cost 96, optimum 79 |
| `data/r3_path_witness_5_servers.json` | Sharp `r=3, K=2` path: 5 clients, 5 servers, cost 6, optimum 5 |
| `data/small_graph_verification.json` | Exhaustive small-graph and tree/path check results |

In the witness JSON files, clients and servers are indexed from zero.
`eligibility[u]` lists the servers available to client `u`;
`current_assignment[u]` and `optimal_assignment[u]` give its assigned servers.

No license has been selected for this repository.
