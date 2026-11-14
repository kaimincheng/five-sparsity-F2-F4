# Exact verification certificates for five-sparsity over F₂ and F₄

This repository contains exact and reproducible calculations supporting the paper
**“The five-sparsity of X^n − 1 over the binary and quaternary fields”**
by **Kaimin Cheng and Cong Mao**.

The paper classifies all positive integers \(n\) such that every irreducible factor
of \(X^n-1\) over \(\mathbb F_2\) or \(\mathbb F_4\) has at most five nonzero
coefficients. Its Section 7 (*Finite verifications and reproducibility*)
registers the finite checks **(C1)–(C9)**. The files here follow the same numbering.

This is a **mathematical proof certificate**, not an experimental statistical
study. All arithmetic uses exact integers and finite fields. In particular,
the rank-two collision classification is exhaustive over all candidate
integer-affine lines and does **not** truncate unbounded exponents.

## How to reproduce

On a machine with **Python 3.8+**, **GNU C++17** (`g++`), and a shell:

```bash
bash run_all.sh
```

The script compiles the bounded-degree enumerator, performs each check,
terminates immediately if an assertion fails, and saves the combined output to
[`logs/verification.log`](logs/verification.log). No third-party Python packages
are needed. On Windows, a terminal with `g++`, Python, and Bash (for instance
WSL or Git Bash with a working C++ compiler) can run the same command.

If a C++17 compiler is named differently, set `CXX` before running, for
example `CXX=clang++ bash run_all.sh`.

## Verification register

| Check | Mathematical assertion verified | Source code / data |
|---|---|---|
| C1 | For odd \(15\le t\le325\), exactly **2,821,078** candidate binary quintinomials; exactly **272** have \(\mathrm{wt}(\mathcal N_3(f))\le5\), and their minimum fifth-norm weight is **7**. | [`verify_R_finite.cpp`](verify_R_finite.cpp) |
| C2 | For the 27 residue triples, formal cubic support sizes **13, 14, 17** occur **18, 7, 2** times; exactly **75** residue-labelled collision hyperplanes have generic weight at most five. | [`verify_collision_lattices.py`](verify_collision_lattices.py), [`verify_collision_lattices_alt.py`](verify_collision_lattices_alt.py) |
| C3 | Exactly **11,298** integral affine lines arise from independent collision hyperplanes; **1,026** have generic weight at most five; only **seven** pass the unrestricted admissibility/parity test. | [`verify_collision_lattices.py`](verify_collision_lattices.py), [`verify_collision_lattices_alt.py`](verify_collision_lattices_alt.py) |
| C4 | The 35 fifth-transform supports of the seven infinite families contain **338** affine exponent pairs; determinant and conjugate-product methods agree, and no further collision occurs for any integer \(K\ge1\). | [`verify_fifth_norm.py`](verify_fifth_norm.py), [`data/norm5_affine_supports.csv`](data/norm5_affine_supports.csv) |
| C5 | Exact order-prime lists for odd orders below 15; irreducible factors of weights 7, 9, 7 at \(p=23,89,127\). | [`verify_prime_exceptions.py`](verify_prime_exceptions.py) |
| C6 | Explicit factorization of \(\Phi_{13}\) into two degree-six factors over \(\mathbb F_4\). | [`verify_prime_exceptions.py`](verify_prime_exceptions.py) |
| C7 | The eight listed binary mixed-order witnesses are irreducible factors of \(\Phi_d\) with more than five terms. | [`verify_binary_mixing.py`](verify_binary_mixing.py) |
| C8 | The 16 quaternary obstruction certificates: prescribed coefficient of the product of elementary-symmetric expressions reduced modulo \(\Phi_d\) is 1. | [`verify_quaternary_mixing.py`](verify_quaternary_mixing.py) |
| C9 | Exact factorization and weights of \(X^d-1\) for \(d=3,5,7,9,15,21,31,63,73\); in particular the eight degree-nine factors of \(\Phi_{73}\) and the two stated factors of \(\Phi_{21}\). | [`verify_binary_cores.py`](verify_binary_cores.py) |

The shared utility [`gf2.py`](gf2.py) implements polynomial arithmetic over
\(\mathbb F_2\) by integer bitmasks, with Rabin's irreducibility test.

## Methods and auditability

- **C1** exhausts the mathematically bounded degree range. The third norm
  uses the cubic resultant identity and the fifth norm uses a product of
  fifth-root conjugates over \(\mathbb F_{16}\).
- **C2–C3** use **two independent implementations**. They construct formal affine exponent functions symbolically,
  normalize collision hyperplanes, intersect independent pairs, test the
  integer-lattice criterion via the cross product and Bézout's identity,
  and solve linear inequalities plus parity over a full integer line.
  The independent checker uses 2-by-2 and 3-by-3 exterior minors,
  a minor-modulus integer-point test, and projective Pluecker coordinates.
  No bound on \(T,A,B\) or the line parameter is introduced.
- **C4** independently computes every fifth norm as (i) the characteristic-two
  determinant of a 5-by-5 multiplication matrix, and (ii) a product of
  conjugates in \(\mathbb F_{16}\). These yield exactly the same complete
  supports; a rational divisibility check rules out all pairwise collisions
  for integer \(K\ge1\).
- **C5, C7, C9** use bitwise arithmetic in \(\mathbb F_2[X]\), exact
  multiplicative-order tests, Rabin's irreducibility criterion, and exact
  product identities.
- **C6** multiplies two polynomials over
  \(\mathbb F_4=\mathbb F_2[\omega]/(\omega^2+\omega+1)\).
- **C8** works entirely in
  \(\mathbb F_2[T]/(T^d-1)\) to construct the elementary-symmetric
  polynomials and reduces their product modulo the binary reduction of
  \(\Phi_d\). The latter is computed recursively by *exact polynomial
  division*, without a computer algebra package.

## Exactness and software

The code uses only Python's standard library and a C++17 compiler. It does
not use numerical approximation, heuristic search, randomization, or
proprietary libraries. Every claimed result is guarded by an assertion in
its corresponding program. The `data/` directory stores the **full** list
of the 338 verified fifth-transform affine exponent pairs; the log stores
an example of a successful complete run.

## Provenance and citation

This repository contains an independently implemented set of checking
programs prepared for the above manuscript. Its numbered outputs match the
finite verification register in Section 7. The repository should be
associated with the submitted manuscript version and pinned by a Git tag or
release; do not change published certificate results without recording a
new version.

When citing this software, include the repository URL and a permanent tag
or commit identifier. This document deliberately does not invent a DOI,
release number, or journal acceptance date.
