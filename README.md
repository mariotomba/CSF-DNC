# Chromatic symmetric functions via deletion-near-contraction

This repository contains SageMath research code for computing the chromatic symmetric function (CSF) in the star basis. The implementation uses the deletion-near-contraction (DNC) relation introduced by Aliste-Prieto, de Mier, Orellana, and Zamora.

These notebooks contain companion research code for the published paper *The chromatic symmetric function in the star-basis* by Michael Gonzalez, Rosa Orellana, and Mario Tomba.

## Requirements

- [SageMath](https://www.sagemath.org/) with a Jupyter kernel
- NumPy (available in standard SageMath installations)

The notebooks were created with SageMath 9.5 and 10.6 kernels. SageMath 10.6 is the most recent version recorded in the repository and is the recommended starting point.

## Quick start

Clone the repository, change into the notebook directory, and start Jupyter through SageMath:

```bash
git clone https://github.com/mariotomba/CSF-DNC.git
cd CSF-DNC/notebooks
sage -n jupyter
```

Open a notebook and run its cells from top to bottom. The computational notebooks load the shared definitions with:

```python
%run DNCfundamentals.ipynb
%run helperFunctions.ipynb
```

Run the notebooks from the `notebooks/` directory so these relative paths resolve correctly.

## Notebook guide

| Notebook | Purpose |
| --- | --- |
| `DNCfundamentals.ipynb` | Core deletion-near-contraction recurrence, CSF computation, graph constructors, and display helpers. |
| `helperFunctions.ipynb` | Shared utilities for partitions, printed polynomials, and caterpillar data. |
| `LeafContractionPolynomial.ipynb` | Computes the leaf-contraction polynomial and extracts its terms from a CSF vector. |
| `CSFmatrix-caterpillar-basis.ipynb` | Builds the caterpillar-basis CSF matrix and computes basis coefficients. |
| `CSFposet.ipynb` | Experimental displays of CSF-supported partitions grouped by length or number of ones. |
| `CaterpillarsExperiments.ipynb` | Exploratory computations for caterpillars and comparisons of CSF-derived invariants. |

## Minimal example

After loading the two shared notebooks, construct a caterpillar from a composition and compute its CSF:

```python
T = create_caterpillar([2, 3, 1, 2, 3, 4])
csf_vector = CSF_tree(T, prnt=1)
```

The returned vector is indexed by `Partitions(n)` in SageMath's default order, where `n` is the number of vertices.

## Trees with equal-order leaf components

In `DNCfundamentals.ipynb`, run the definition cells through **Equal-order leaf
components**, then call this function in a new cell:

```python
CSF_uniform_leaf_components(8, 2)
```

The two inputs are the total number of vertices `n` and the common leaf-component
order `l`, both integers at least 2. A leaf component is a connected component
obtained by deleting all internal edges (edges whose endpoints both have degree
greater than one).

For each tree, once up to isomorphism, the function displays its plot, prints
`Leading partition: ...`, and prints the CSF as a linear combination of `st[...]`,
using the same DNC computation and formatting as `CSF_iterator`. The example above
displays two trees, both with leading partition `[2, 2, 2, 2]`.

There are no qualifying trees if `l` does not divide `n`; the function prints a
message in that case. Otherwise it generates cores on `n/l` vertices and attaches
`l-1` leaves to every core vertex, avoiding enumeration of all trees on `n`
vertices. When `n == l`, it displays the single star, including `K2` for `(2, 2)`.
The number of cores and the DNC computation can still grow quickly.

To use the trees and coefficient vectors without displaying them:

```python
for T in trees_uniform_leaf_components(12, 3):
    csf_vector = CSF_tree(T)
    # Process T and csf_vector here.
```

Run the regression tests from the repository root with:

```bash
sage -python -m unittest discover -s tests -v
```

## Reproducibility notes

- Some exploratory cells enumerate all trees or many caterpillar compositions. Their runtime and memory use grow quickly with `n`; begin with smaller values when checking an installation.
- Notebook outputs have been retained because they record examples used during the research process.
- `DNCfundamentals.ipynb` and `helperFunctions.ipynb` define the shared functions used by the remaining notebooks.

## Publication and citation

This repository accompanies the following published paper:

> Michael Gonzalez, Rosa Orellana, and Mario Tomba. “The chromatic symmetric function in the star-basis.” *Discrete Mathematics* 349, no. 2 (2026), Article 114691. [https://doi.org/10.1016/j.disc.2025.114691](https://doi.org/10.1016/j.disc.2025.114691)

An author manuscript is also available on [arXiv:2404.06002](https://arxiv.org/abs/2404.06002).

If you use this code, please cite the published paper above as well as this repository.

## Acknowledgments

This code was developed as part of joint work with Rosa Orellana and Michael Gonzalez.

The DNC recurrence implemented here is due to Aliste-Prieto, de Mier, Orellana, and Zamora. The caterpillar-basis and related experiments belong to the Gonzalez–Orellana–Tomba project described above.
