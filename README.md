# symmetric-subgraph-detection

Header-only C++17 implementation of a shifted-LDL^T symmetry detector for
undirected graphs. Given a graph and a candidate shift `mu`, the detector returns
vertex groups that support vectors of the adjacency eigenspace near `mu`.

The method comes from Algorithm 1 of *Formulas for Detecting Symmetries in Graphs
via Shifted Factorization*. For a disjoint union of cliques of size at least two
and `mu = -1`, the groups are exactly the cliques. For a 4-cycle and `mu = 0`,
this implementation returns the two opposite-vertex pairs; the partition depends
on the basis, as described in [docs/algorithm.md](docs/algorithm.md). The
detector is not an automorphism enumerator; see the same note for the precise
statement and its limits.

The project was originally a Python notebook built on NetworkX and SciPy. This
repository is a dependency-free C++ rewrite with its own LDL^T factorization, a
test suite, and a command-line driver.

## Build and test

Requirements: a C++17 compiler and CMake 3.16 or newer.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd build && ctest --output-on-failure
```

The library itself has no build step. Add `include/` to the include path, or add
this directory with `add_subdirectory()` and link `sna::sna`. No installed CMake
package config is provided.

## Quick start

```cpp
#include <iostream>
#include <sna/sna.hpp>

int main() {
    sna::Graph<double> graph;
    for (int i = 0; i < 10; ++i) {
        for (int j = i + 1; j < 10; ++j) graph.addEdge(i, j);
    }

    for (const auto& group : sna::detectSymmetries(graph)) {
        for (int vertex : group) std::cout << vertex << " ";
        std::cout << "\n";
    }
}
```

For the notebook graph (five disjoint 10-cliques plus isolated vertices),
`examples/cliques.cpp` prints:

```
vertices: 53, edges: 225
groups: 5
  10 vertices: 0 1 2 3 4 5 6 7 8 9
  10 vertices: 10 11 12 13 14 15 16 17 18 19
  10 vertices: 20 21 22 23 24 25 26 27 28 29
  10 vertices: 30 31 32 33 34 35 36 37 38 39
  10 vertices: 40 41 42 43 44 45 46 47 48 49
LDL^T relative residual: 0
```

## What the detector does

For `B = A - mu*I`, the code factors `B` as `P B P^T = L*D*L^T` with 1x1 and
2x2 pivots, and reports `L`, the blocks of `D`, and the permutation. A near-zero 1x1 pivot, or a near-zero eigenvalue of a 2x2 pivot,
yields a seed; solving `L^T v = seed` and undoing the permutation gives an
approximate null vector. The supports of these vectors are joined and split into
connected components. Groups with fewer than two vertices are dropped.

The complete-pivoting rule in `factorSymmetric` guarantees that a zero pivot is
emitted only when the remaining block is zero, so the exact-arithmetic
factorization never divides by zero even though `B` is singular by construction.
The derivation is in [docs/algorithm.md](docs/algorithm.md).

## API

All functions are in `namespace sna` and are available through `sna/sna.hpp`.

### Detection

| Function | Result |
| --- | --- |
| `detectSymmetries(const Matrix<T>&, const DetectorParams<T>& = {})` | groups of row/column indices |
| `detectSymmetries(const Matrix<T>&, const std::vector<T>& shifts, ...)` | union of per-shift groups, deduplicated |
| `detectSymmetries(const Graph<T>&, const DetectorParams<T>& = {})` | groups of vertex labels |
| `detectSymmetries(const Graph<T>&, const std::vector<T>& shifts, ...)` | union of per-shift groups, deduplicated |

`DetectorParams<T>` has three fields:

```cpp
T shift       = -1;     // mu
T block_tol   = 1e-8;   // paper epsilon: near-zero block eigenvalue
T support_tol = 1e-5;   // paper tau: support threshold on a unit null vector
```

Groups contain at least two vertices. Vertices within a group are sorted
ascending and the list of groups is sorted lexicographically, so output is
deterministic. A matrix overload returns indices because a bare matrix has no
vertex labels; a graph overload returns the labels passed to `Graph::addEdge` or
`Graph::addNode`.

The matrix overloads assume symmetry and check only that the matrix is square.
Call `isSymmetric()` if the input is not under your control. For a multi-shift
call, the `shift` field of the supplied params is overwritten for each shift; the
result can contain overlapping groups because it is a union, not a partition.

### Numerical pieces

| Function | Purpose |
| --- | --- |
| `factorSymmetric(const Matrix<T>&)` | `L`, the 1x1/2x2 blocks of `D` (`PivotBlock`), and `perm` such that `A(perm[i],perm[j]) = (L D L^T)(i,j)` |
| `reconstructionResidual(A, factor)` | relative Frobenius residual of `A(perm[i],perm[j]) = (L D L^T)(i,j)` |
| `solveTransposedUnitLower(L, rhs)` | solves `L^T x = rhs` for unit lower triangular `L` |
| `approxNullSpace(A, params)` | unit-norm basis of the approximate null space of `A - shift*I` |
| `groupBySupport(vectors, tol)` | support union and connected components |
| `nullSpaceResidual(A, shift, x)` | `‖(A - shift*I)x‖ / ‖x‖`; returns 0 for `x = 0` |
| `inducedDensity(A, vertices)` | weighted density of an induced submatrix |

`approxNullSpace` returns unit vectors, not an orthogonal basis.
`groupBySupport` requires all vectors to have the same length and uses them as
given: it does not normalize, so pass unit vectors if `support_tol` is meant as a
relative threshold. It drops singleton supports.

### Refinement

| Function | Result |
| --- | --- |
| `filterByDensity(A or Graph, groups, min_density)` | keeps groups whose weighted density is at least `min_density` |

`filterByDensity` uses `inducedDensity`, removes duplicate vertices before
computing the density, and returns the cleaned groups. It filters; it does not
split a group or resolve overlapping supports.

### Containers

`Matrix<T>` is a dense row-major matrix; element access is unchecked. `Graph<T>`
is an undirected weighted graph with integer labels. `Graph::addEdge` creates
missing endpoints, sums parallel edges, and ignores self-loops (the label is
still created). Edges with weight 0 are stored and appear in `neighbors()`, but
they contribute nothing to `adjacencyMatrix()` or `density()`.
`Graph::adjacencyMatrix()` uses labels sorted ascending, so the result does not
depend on insertion order.

## Tolerances and scale

`block_tol` is absolute. For a matrix whose largest entry in absolute value is
`s`, use about `1e-8 * s`; the default `1e-8` is meant for 0/1 adjacency
matrices. `support_tol` is applied to unit vectors: lowering it includes more
vertices, raising it can split a group.

`mu` must be near an eigenvalue of interest, usually a repeated one. Common
choices:

| Graph | Repeated eigenvalue |
| --- | --- |
| `K_n` | `-1` |
| `K_n` with uniform weight `w` | `-w` |
| `K_{a,b}` | `0` |
| `C_4` | `0` |

## Limits

* The implementation is dense: O(n^3) time and O(n^2) memory. It is not suitable
  for very large sparse graphs, despite the sparse-factorization motivation of
  the original note.
* The detector finds support groups of an eigenspace, not automorphism orbits.
  A clique embedded in a larger graph keeps its full `-1` eigenspace only if the
  outside connections are constant on the clique; otherwise a subspace can
  survive (see `testPartialExternalCoupling`). `testAttachedClique` is the
  uniform case.
* Overlapping symmetric subgraphs merge when their supports connect. No
  max-clique or community split is attempted; the split is ambiguous without a
  model of which symmetry is wanted.
* A simple (non-repeated) eigenvalue also produces a group. For a delocalized
  eigenvector the group can be large and may not look like a motif.
* A vertex that appears only in a one-vertex support is not returned. At
  `mu = 0`, isolated vertices are in the null space but are dropped for this
  reason.

## Performance

`bench/bench.cpp` measures `detectSymmetries` (factorization plus grouping) on a
disjoint union of 10-cliques, one per ten vertices, with the default shift.
Median of five runs in a Release build on an Apple-silicon laptop:

| Vertices | Seconds | Groups |
| --- | --- | --- |
| 100 | 0.0006 | 10 |
| 200 | 0.0045 | 20 |
| 400 | 0.036 | 40 |
| 800 | 0.32 | 80 |

The growth is roughly cubic over this range, as expected from the dense
factorization.

## Repository layout

```
include/sna/      headers: matrix, ldlt, graph, spectral, umbrella
examples/         cliques, extensions (weighted/multi-shift/filter), edge-list CLI
tests/            test binary and helpers
bench/            scaling check
docs/algorithm.md derivation, proofs, traceability, limits
```

Build targets: `example_cliques`, `example_extensions`, `example_edge_list`,
`sna_tests`, `bench_scale`.

The edge-list driver accepts `u v` or `u v w` per line and `#` comments:

```sh
./build/example_edge_list graph.edges -1
```

Vertices that appear in no edge cannot be declared in this format and are not
added to the graph. A non-numeric weight is rejected.

## Tests

`sna_tests` runs checks for:

* LDL^T reconstruction on random symmetric and rank-deficient matrices, with
  worst relative residuals of `3.0e-16` and `5.6e-16`;
* the 2x2 pivot branch, block coverage, unit lower triangularity, permutation
  validity, determinism, and singular pivot cases;
* null-vector residuals and the expected partition on random disjoint clique
  unions, including relabeled graphs;
* a near-zero 2x2 block with one and with two eigenvalues below `block_tol`;
* weighted cliques, a 4-cycle at zero, an attached clique, partial external
  coupling, multi-shift deduplication, the density filter,
  `float`/`double`/`long double`, and degenerate inputs;
* rejection of malformed inputs by the public helpers.

The current test binary reports 18,495 checks with zero failures.

## History

The initial prototype was a Python notebook that used `networkx` and
`scipy.linalg.ldl`. The C++ code keeps the same shift-and-factorization idea but
implements the factorization itself with complete pivoting. The reason is not
that partial pivoting is wrong for nonsingular matrices; it is that this
factorization must remain complete when `B` is singular, and complete pivoting
gives a simple guarantee: an exactly zero pivot is emitted only when the
remaining block is zero.

The prototype and the note are not redistributed in this repository.

## License

MIT. See [LICENSE](LICENSE).
