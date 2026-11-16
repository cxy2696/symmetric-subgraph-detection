# Algorithm and implementation notes

These notes state what the detector computes, why the LDL^T factorization
produces null vectors, and which parts of the method are exact. They are the
specification implemented by this package. The revised working note in
[Finding_symmetries_in_networks.pdf](Finding_symmetries_in_networks.pdf)
motivates the shifted-factorization approach; it is a working note with open
questions, not a fixed algorithm. `block_tol` is the near-zero tolerance and
`support_tol` is the support threshold. The tolerance is relative to `max|B|`,
so the numerical-zero test does not depend on the scale of the matrix.

## 1. Output

Let `A` be the adjacency matrix of an undirected graph, let `mu` be the shift,
and let `B = A - mu*I`. Let `X = ker(B)` be the null space (in exact arithmetic).

For a vector `x`, define its support `supp(x) = { i : |x_i| > tau }`. Pick a basis
`x_1, ..., x_q` of `X`, put a clique on each `supp(x_j)`, and take the connected
components of the resulting support graph. The detector returns those components
that contain at least two vertices. In these notes, "group" means such a component.

The groups are vertex sets that support vectors of the `mu`-eigenspace. They are
not automorphism orbits, and the method does not enumerate symmetries. For a
disjoint union of cliques of size at least two at `mu = -1`, the groups are
exactly the cliques; isolated vertices (`K_1`) are dropped by the two-vertex
rule. That is the case the original notebook exercised.

## 2. From a near-zero block of D to a null vector

The factorization is

```
P B P^T = L D L^T
```

where `P` is a permutation, `L` is unit lower triangular, and `D` is block
diagonal with 1x1 and 2x2 blocks. Write `B_p = L D L^T` for the permuted matrix.

**Claim 1.** Let `u` be a unit vector supported on one diagonal block of `D`,
chosen so that `D u = 0`. Then `v = L^{-T} u` satisfies `B_p v = 0`, and undoing
the permutation gives a null vector of `B`.

*Proof.* `L^T` is invertible, so the system `L^T v = u` has the unique solution
`v = L^{-T} u`. Then `B_p v = L D L^T v = L D u = 0`. If `P` is the permutation
matrix with `B_p = P B P^T`, then `B (P^T v) = 0`. In code, `v` is indexed by
permuted position and `perm[i]` maps position `i` back to the original vertex.

For a 1x1 zero block at position `k`, `u = e_k`. For a 2x2 block `E` with an
eigenvalue `lambda = 0`, `u` is the corresponding eigenvector of `E`, placed at
positions `k` and `k+1`. The 2x2 case needs the eigenvector: using `e_k` gives
`D e_k = (E e_k, 0, ...)`, which is not zero unless `E` has a zero first column.
The code uses this eigenvector.

**Claim 2.** The vectors obtained from all zero blocks form a basis of `ker(B)`.

*Proof.* The kernel of `D` is spanned by the seeds (zero 1x1 blocks and null
eigenvectors of singular 2x2 blocks). Multiplication by `L^{-T}` is an
isomorphism, so the images form a basis of `ker(B_p)`; applying `P^T` gives a
basis of `ker(B)`.

With the pivot rule of Section 3, an exactly singular `2x2` block does not occur:
every 2x2 pivot has `|det| >= 0.75*M^2 > 0`. Exact nullity therefore appears as
zero 1x1 blocks. The 2x2 seed branch is for *numerically* near-zero blocks: when
an eigenvalue of a 2x2 pivot is below `block_tol`, the corresponding seed gives a
vector that is a null vector up to the size of that eigenvalue.

**Claim 3.** The union of the supports is independent of the basis; the partition
into connected components is not.

*Proof.* A vertex `i` is absent from every vector in `X` exactly when row `i` of
the orthogonal projector onto `X` is zero. That condition depends only on `X`.
The partition, on the other hand, is built from the supports of the particular
basis at hand: two bases can give the same union but different components. For
example, a 4-cycle at `mu = 0` has a null space spanned by vectors on `{0,2}` and
`{1,3}`, but the vector `(1,1,-1,-1)` is also a null vector supported on all
four vertices. The implementation returns the partition induced by the basis
that comes out of the factorization, which for this example is `{0,2}`, `{1,3}`.

For a near-zero block, the construction gives an approximate null vector. If a
1x1 block has value `d`, then `B_p v = d L e_k`; if a 2x2 block has eigenvalue
`lambda`, then `B_p v = lambda L u`. In both cases the residual is the block
value times a factor that depends on `L`. The effective threshold is
`block_tol * max|B|`; it bounds the block value, not `cond(L)`, which is not
checked.

## 3. Pivot strategy

`B` is singular by construction, so the factorization has to be able to emit a
zero pivot. That is safe only when the remaining row and column are zero;
otherwise elimination would divide by zero. Complete pivoting over the trailing
block (the Bunch-Parlett strategy), with a preference for diagonal pivots, gives
this invariant directly. Let

```
M = max |S[i][j]| over i != j,   G = max |S[i][i]|
```

in the trailing block `S`.

* If `G == 0` and `M == 0`, the block is zero; emit zero 1x1 pivots.
* If `G >= 0.5*M`, use a 1x1 pivot on the largest diagonal entry.
* Otherwise use the 2x2 block that contains the largest off-diagonal entry.

The rule has three consequences, in exact arithmetic.

1. If `M = 0`, all off-diagonals are zero and no elimination is needed. If
   `M > 0`, a 1x1 pivot used for elimination is nonzero: `G >= 0.5*M > 0`, and
   every multiplier `|S[i][k] / d|` is bounded by `M / (0.5*M) = 2`.
2. A 2x2 pivot is invertible. Its off-diagonal entry is `b` with `|b| = M`, and
   its diagonal entries satisfy `|a|, |c| < 0.5*M`, so
   `|det| = |a*c - b^2| >= M^2 - (0.5*M)^2 = 0.75*M^2`. The multipliers are
   bounded by 2 as well.
3. A zero pivot occurs only when every entry of the trailing block is zero.
   There is then nothing to eliminate, so a zero 1x1 pivot is legitimate.

The exact-arithmetic factorization therefore never divides by an exactly zero
pivot. Floating-point behavior still depends on the conditioning of `B`; the
bounds above are not a backward-error analysis. The O((n-k)^2) maximum search at
step `k` keeps the total cost at O(n^3). The tests in `tests/test_ldlt.cpp`
include singular matrices and rank-deficient random matrices.

A consequence for the detector: complete pivoting takes large entries first. A
well-scaled 2x2 pivot has `|det| >= 0.75*M^2` and both eigenvalues are
`Omega(M)`. The detector's 2x2 seed branch therefore fires only when the pivot
block's eigenvalue is below `block_tol`, which for a 2x2 pivot means the block's
scale `M` is also `O(block_tol)`. This is the numerical-null-space case. Small
eigenvalues of the original matrix appear in the trailing blocks of the
factorization; they may be 1x1 or 2x2 pivots.

## 4. Steps and cost

| Step | Where | Cost |
| --- | --- | --- |
| Build the adjacency matrix | `Graph::adjacencyMatrix` | O(n^2 + E) |
| Shift by `mu` | `approxNullSpace` | O(n^2) |
| LDL^T | `factorSymmetric` | O(n^3) |
| Scan `D` for near-zero blocks | `approxNullSpace` | O(n) |
| One triangular solve per seed | `solveTransposedUnitLower` | O(q*n^2), q <= n |
| Undo the permutation, normalize | `approxNullSpace` | O(q*n) |
| Support union and components | `groupBySupport` | O(q*n) |
| Density filter | `filterByDensity` | O(sum s^2), s = group size |

Total time is O(n^3), dominated by the factorization. Memory is O(n^2) for the
matrices and the basis. The implementation is dense and targets small and medium
matrices. The asymptotic cost is the same as a dense symmetric eigendecomposition;
this code does not accumulate eigenvectors, but no constant-factor comparison is
measured here.

## 5. Implementation map

| Step | Implementation |
| --- | --- |
| Build `B = A - mu*I` and scale it | `approxNullSpace` |
| Factor `P B P^T = L D L^T` | `factorSymmetric` |
| Scan near-zero blocks of `D` | block loop in `approxNullSpace` |
| Solve `L^T v = u`, normalize, undo `P` | `addSeed` lambda |
| Union supports, take components | `groupBySupport` |
| Density filter | `filterByDensity` |
| Iterate over shifts | `detectSymmetries(matrix, shifts)` |
| Return sorted groups | `detectSymmetries` |

Implementation choices worth recording:

* The support graph places a clique on each support. Connecting every support
  vertex to the first one spans the same components in O(|supp|) instead of
  O(|supp|^2).
* Refinement is not part of `detectSymmetries`. The caller applies a policy, or
  none. The package provides a density filter only; it does not compute a
  maximal clique or split an overlapping support, because that split is not
  determined by the null-space data alone.
* Pivot selection is complete pivoting (Section 3).
* `block_tol` is relative to `max|B|` rather than absolute (Section 7). The
  shifted matrix is divided by `max|B|` before factorization, which leaves the
  null space unchanged and keeps the pivots in a comparable range. For 0/1
  adjacency matrices the two conventions agree in scale.
* Results are sorted deterministically: vertices ascending within a group,
  groups lexicographically.
* Matrix overloads use row/column indices; graph overloads translate back to
  vertex labels.

## 6. Exactness and limits

1. **Numerical, not symbolic.** A block is treated as zero when it is below
   `block_tol * max|B|`, and a vector entry joins a support when it exceeds
   `support_tol` on a unit vector. The first threshold is relative to the
   matrix, the second to the vector; neither depends on the ambient scale. The
   matrix scale is global to `B`: if `B` has blocks with very different scales,
   a near-nullity in a small block can fall below the threshold set by the
   largest block. Run the blocks separately or rescale them in that case.

2. **External coupling.** Let `S` be a unit-weight clique in a larger graph with
   no self-loops, and let `x` be supported on `S` with `sum(x) = 0`. For `i` in
   `S`, `(A x)_i = -x_i`, but for `k` outside `S`,
   `(A x)_k = sum_{i in S} A[k][i] x_i`. The full null space on `S` survives
   only when each outside row is constant on `S`. Otherwise a smaller subspace
   can survive: for a `K4` on `0..3` with an outside vertex joined to `0,1` with
   weight 1 and to `2,3` with weight 2, the vectors on `{0,1}` and `{2,3}` are
   still null vectors. The test `testAttachedClique` is the uniform case: a `K4`
   with an extra vertex joined to all four vertices and a pendant on that extra
   vertex.

3. **Overlapping supports merge.** If two symmetric subgraphs share vertices,
   their supports may connect in the support graph and be returned as one group.
   The package does not split such a group; `filterByDensity` can reject it, but
   separating it requires a model of which symmetry is wanted.

4. **No multiplicity requirement.** A simple eigenvalue also produces a null
   vector and therefore a group. The method does not test whether the eigenvalue
   is repeated. For a simple eigenvalue the vector is often delocalized, and the
   returned group can be large and not correspond to a visible symmetry.

5. **Singleton supports.** A vertex that appears only in a one-vertex support is
   not returned. This makes the returned union basis-dependent even though the
   active set of `X` is not. The zero matrix at `mu = 0` is an extreme example:
   the basis from the factorization is the standard basis, every support is a
   singleton, and the detector returns no groups, although a basis such as
   `e_0 + e_1` would produce a group.

6. **Dense cost.** O(n^3) time and O(n^2) memory. The benchmark in `bench/`
   shows the expected growth. Graphs with tens of thousands of vertices need a
   sparse factorization, which is outside this package.

## 7. Choosing `mu`, `block_tol`, and `support_tol`

`mu` must be near an eigenvalue of interest, usually a repeated one. Useful
cases:

| Graph | Spectrum | Repeated value |
| --- | --- | --- |
| `K_n`, unit weights | `n-1`, `-1` (n-1 times) | `-1` |
| `K_n`, uniform weight `w` | `w(n-1)`, `-w` (n-1 times) | `-w` |
| `K_{a,b}` | `+sqrt(ab)`, `-sqrt(ab)`, `0` (a+b-2 times) | `0` |
| `C_4` | `2`, `0`, `0`, `-2` | `0` |

`block_tol` is relative to `max|B|`: a pivot or 2x2 eigenvalue is zero when its
magnitude is below `block_tol * max|B|`. The default `1e-8` is intended for 0/1
adjacency matrices, where `max|B|` is of order one for a shift of order one.
This margin is far above the rounding of an exact zero pivot, which scales with
machine epsilon; it is a tolerance, not a formal error bound.

`support_tol` is applied to unit vectors. The default `1e-5` separates genuine
components from rounding noise for the graphs in the tests. Lowering it admits
more vertices; raising it can split a genuine group.

## 8. Validation

`tests/` checks the following:

* LDL^T reconstruction on random symmetric matrices and on rank-deficient
  matrices, with worst relative residuals of `3.0e-16` and `5.6e-16`.
* The 2x2 pivot branch, unit lower triangularity, block coverage, permutation
  validity, determinism, and singular pivot cases.
* Null-vector residuals: the two exactly singular test matrices (`K_5` and the
  3-vertex path) check that every returned vector has residual at most `1e-12`;
  the random clique-union tests check the expected group partition.
* A near-zero 2x2 block where one eigenvalue, and separately both eigenvalues,
  fall below `block_tol * max|B|`; this exercises the 2x2 seed branch.
* The same structure at scales `1e-20`, `1` and `1e6`, which checks that the
  relative tolerance is scale-invariant; non-finite entries and negative
  tolerances are rejected.
* The notebook example (`K_{10}` five times plus isolated vertices), weighted
  cliques, a 4-cycle at zero, an attached clique, multi-shift deduplication,
  the density filter, scalar types (`float`, `double`, `long double`), and
  empty/degenerate inputs. Public helpers reject malformed input instead of
  reading out of bounds.
