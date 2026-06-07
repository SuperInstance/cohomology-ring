# cohomology-ring

> **Cup product. Cohomology operations. The ring structure of H*(X).**

[![crates.io](https://img.shields.io/crates/v/cohomology-ring.svg)](https://crates.io/crates/cohomology-ring)
[![docs.rs](https://docs.rs/cohomology-ring/badge.svg)](https://docs.rs/cohomology-ring)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust library for computing cohomology rings and cohomology operations. Implements cochain complexes, cup product (giving cohomology its multiplicative ring structure), Bockstein homomorphism, and Steenrod squares. Distinguishes spaces that homology alone cannot — the algebraic structure that makes cohomology strictly more powerful than homology.

---

## Table of Contents

- [What is a Cohomology Ring?](#what-is-a-cohomology-ring)
- [Why Does This Matter?](#why-does-this-matter)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
- [Mathematical Background](#mathematical-background)
- [Installation](#installation)
- [Related Crates](#related-crates)
- [License](#license)

---

## What is a Cohomology Ring?

**Cohomology** assigns to each topological space X a sequence of groups H⁰(X), H¹(X), H²(X), ... that capture the space's structure. But cohomology is more than groups — it's a **graded ring**.

The **cup product** ∪: Hᵏ × Hˡ → Hᵏ⁺ˡ gives cohomology a multiplicative structure:

```
α ∈ Hᵏ(X),  β ∈ Hˡ(X)  →  α ∪ β ∈ Hᵏ⁺ˡ(X)
```

This ring structure is a strictly stronger invariant than cohomology groups alone. For example:
- CP² (complex projective plane) and S² ∨ S⁴ (wedge of spheres) have identical cohomology **groups**
- But their cohomology **rings** differ: in CP², α ∪ α ≠ 0 for α ∈ H², while in S² ∨ S⁴, α ∪ α = 0

**Cohomology operations** go further:
- **Bockstein** β: Hᵏ(X; Z/p) → Hᵏ⁺¹(X; Z/p) — detects torsion
- **Steenrod squares** Sqⁱ: Hᵏ(X; Z/2) → Hᵏ⁺ⁱ(X; Z/2) — stable cohomology operations

These operations provide information that even the ring structure alone cannot capture.

## Why Does This Matter?

**For topology**: Cohomology rings are fundamental invariants for classifying spaces. Two spaces with the same cohomology ring are "close" topologically — the ring is a powerful discriminator.

**For data analysis**: Applied topology (TDA) uses cohomology to analyze the shape of data. The cup product captures relationships between features in different dimensions — how 1D loops combine to create 2D voids.

**For physics**: Cohomology rings classify characteristic classes of vector bundles — the mathematical language of gauge theories, magnetic monopoles, and topological insulators.

**For agent systems**: Topological invariants of agent state spaces classify the "shape" of agent behavior. The ring structure captures how behavioral primitives compose into complex patterns.

## Architecture

```
cohomology-ring
│
├── Cochain                    ← Basic cochain element
│   ├── zero(degree)               Zero cochain
│   ├── basis(degree, index, coeff) Single basis element
│   ├── add(&other)                Cochain addition
│   ├── scale(scalar)              Scalar multiplication
│   ├── is_zero()                  Triviality check
│   └── entries()                  Sparse coefficient access
│
├── CochainComplex             ← Sequence of cochain groups + coboundary
│   ├── new(dimensions, coboundary) Custom complex
│   ├── trivial(dimensions)        Zero coboundary
│   ├── for_circle()               S¹ complex
│   ├── for_sphere()               S^n complex
│   ├── for_rp2()                  RP² complex
│   ├── for_torus()                T² complex
│   ├── apply_coboundary(cochain)  δ: Cᵏ → Cᵏ⁺¹
│   └── coboundary_matrix(k)       Matrix representation of δ
│
├── CohomologyGroup            ← Computed cohomology
│   ├── compute(complex, k)        Hᵏ = ker(δᵏ)/im(δᵏ⁻¹)
│   ├── compute_all(complex)       All H⁰, H¹, ...
│   ├── dimension                  Rank of cohomology group
│   └── generators                 Basis representatives
│
├── CupProduct                 ← Ring multiplication
│   ├── cup(a, b)                  α ∪ β
│   ├── verify_graded_commutativity()  α ∪ β = (−1)^kl β ∪ α
│   └── verify_associativity()     (α ∪ β) ∪ γ = α ∪ (β ∪ γ)
│
├── BocksteinHomomorphism      ← Torsion detection
│   ├── apply(cochain, complex)    β: Hᵏ(Z/p) → Hᵏ⁺¹(Z/p)
│   └── verify_square_zero()       β² = 0
│
└── SteenrodSquare             ← Stable operations (mod 2)
    ├── sq0(cochain)               Sq⁰ = identity
    ├── sqn(cochain)               Sqⁿ for top class
    ├── sq1(cochain, complex)      Sq¹ = Bockstein (mod 2)
    ├── verify_cartan_formula()    Sqᵏ(α∪β) = Σ Sqⁱ(α)∪Sqᵏ⁻ⁱ(β)
    └── adem_relation_check(i, j)  Adem relations
```

## Quick Start

```rust
use cohomology_ring::{
    Cochain, CochainComplex, CupProduct,
    BocksteinHomomorphism, SteenrodSquare,
};

// Use a pre-built complex for the torus T²
let torus = CochainComplex::for_torus();

// Compute all cohomology groups
let groups = CohomologyGroup::compute_all(&torus);
for (k, group) in groups.iter().enumerate() {
    println!("H^{}(T²): dimension = {}", k, group.dimension);
}
// H⁰ = 1, H¹ = 2, H² = 1, H³ = 0

// Compute cup product
let a = Cochain::basis(1, 0, 1);  // α ∈ H¹
let b = Cochain::basis(1, 1, 1);  // β ∈ H¹
let ab = CupProduct::cup(&a, &b);  // α ∪ β ∈ H²
println!("α ∪ β: {:?}", ab);

// Verify graded commutativity: α ∪ β = -(β ∪ α) for odd-degree cochains
let is_commutative = CupProduct::verify_graded_commutativity(&a, &b);
println!("Graded commutative: {}", is_commutative);

// Bockstein homomorphism
let x = Cochain::basis(1, 0, 1);
let beta_x = BocksteinHomomorphism::apply(&x, &torus);
println!("β(x): {:?}", beta_x);

// Verify β² = 0
let is_nilpotent = BocksteinHomomorphism::verify_square_zero(&x, &torus);
println!("β² = 0: {}", is_nilpotent);

// Steenrod squares
let y = Cochain::basis(1, 0, 1);
let sq1_y = SteenrodSquare::sq1(&y, &torus);
println!("Sq¹(y): {:?}", sq1_y);

// Verify Cartan formula: Sqᵏ(α∪β) = Σᵢ Sqⁱ(α) ∪ Sqᵏ⁻ⁱ(β)
let cartan_ok = SteenrodSquare::verify_cartan_formula(&a, &b, &torus);
println!("Cartan formula holds: {}", cartan_ok);
```

## API Reference

### Cochain

| Method | Returns | Description |
|--------|---------|-------------|
| `zero(degree)` | `Self` | Zero cochain |
| `basis(degree, index, coeff)` | `Self` | Single basis element |
| `add(&other)` | `Cochain` | Cochain addition |
| `scale(scalar)` | `Cochain` | Scalar multiplication |
| `is_zero()` | `bool` | Triviality check |
| `get(index)` | `i32` | Coefficient at index |
| `entries()` | `&HashMap<usize, i32>` | All non-zero entries |

### CochainComplex

| Method | Returns | Description |
|--------|---------|-------------|
| `new(dims, coboundary)` | `Self` | Custom complex |
| `trivial(dims)` | `Self` | Zero coboundary |
| `for_circle()` | `Self` | S¹ |
| `for_sphere()` | `Self` | S^n |
| `for_rp2()` | `Self` | Real projective plane |
| `for_torus()` | `Self` | T² |
| `apply_coboundary(cochain)` | `Cochain` | δ(cochain) |
| `coboundary_matrix(k)` | `Vec<Vec<i32>>` | Matrix for δᵏ |

### CohomologyGroup

| Method | Returns | Description |
|--------|---------|-------------|
| `compute(complex, k)` | `Self` | Hᵏ = ker δᵏ / im δᵏ⁻¹ |
| `compute_all(complex)` | `Vec<Self>` | All groups H⁰, H¹, ... |

### CupProduct

| Method | Returns | Description |
|--------|---------|-------------|
| `cup(a, b)` | `Cochain` | α ∪ β ∈ Hᵏ⁺ˡ |
| `verify_graded_commutativity(a, b)` | `bool` | α∪β = (−1)ᵏˡ β∪α |
| `verify_associativity(a, b, c)` | `bool` | (α∪β)∪γ = α∪(β∪γ) |

### Bockstein & Steenrod

| Method | Returns | Description |
|--------|---------|-------------|
| `BocksteinHomomorphism::apply(x, complex)` | `Cochain` | β(x) |
| `BocksteinHomomorphism::verify_square_zero(x, complex)` | `bool` | β² = 0 |
| `SteenrodSquare::sq0(x)` | `Cochain` | Sq⁰ = id |
| `SteenrodSquare::sq1(x, complex)` | `Cochain` | Sq¹ |
| `SteenrodSquare::verify_cartan_formula(x, y, complex)` | `bool` | Cartan formula |
| `SteenrodSquare::adem_relation_check(i, j)` | `bool` | Adem relations |

## Mathematical Background

### Cohomology

Given a cochain complex (C*, δ) with δᵏ: Cᵏ → Cᵏ⁺¹ and δ² = 0:

```
0 → C⁰ →δ₀ C¹ →δ₁ C² →δ₂ C³ → ...
```

The k-th cohomology group: Hᵏ = ker(δᵏ) / im(δᵏ⁻¹) = Zᵏ/Bᵏ

### Cup Product

The cup product ∪: Cᵏ × Cˡ → Cᵏ⁺ˡ is defined on simplicial cochains by:

```
(α ∪ β)([v₀, ..., v_{k+l}]) = α([v₀, ..., vₖ]) · β([vₖ, ..., v_{k+l}])
```

Properties:
- **Graded commutativity**: α ∪ β = (−1)ᵏˡ (β ∪ α)
- **Associativity**: (α ∪ β) ∪ γ = α ∪ (β ∪ γ)
- **Distributivity**: α ∪ (β + γ) = α ∪ β + α ∪ γ

### Bockstein Homomorphism

For the short exact sequence 0 → Z → Z → Z/p → 0, the Bockstein β is the connecting homomorphism:

```
β: Hᵏ(X; Z/p) → Hᵏ⁺¹(X; Z/p)
```

Key property: **β² = 0** (it's a differential). The Bockstein detects p-torsion in integral cohomology.

### Steenrod Squares

Sqⁱ: Hᵏ(X; Z/2) → Hᵏ⁺ⁱ(X; Z/2) are stable cohomology operations satisfying:

- **Sq⁰** = identity
- **Sqᵏ(x)** = x² for x ∈ Hᵏ (for x of degree k)
- **Cartan formula**: Sqᵏ(x ∪ y) = Σᵢ₌₀ᵏ Sqⁱ(x) ∪ Sqᵏ⁻ⁱ(y)
- **Adem relations**: SqⁱSqʲ = Σₖ C(j−k−1, 2k−i) Sqⁱ⁺ʲ⁻ᵏSqᵏ for i < 2j

The Steenrod algebra (generated by all Sqⁱ) is a fundamental tool in homotopy theory and stable homotopy.

## Installation

```bash
cargo add cohomology-ring
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
cohomology-ring = "0.1"
```

## Related Crates

Part of the **SuperInstance Exocortex** math fleet:

- **[graph-homology](https://github.com/SuperInstance/graph-homology)** — Clique complexes and Betti numbers
- **[sheaf-laplacian](https://github.com/SuperInstance/sheaf-laplacian)** — Sheaf Laplacian and Hodge decomposition
- **[persistent-agent](https://github.com/SuperInstance/persistent-agent)** — Persistent homology for agent behavior
- **[tropical-graph](https://github.com/SuperInstance/tropical-graph)** — Max-plus algebra on graphs
- **[categorical-coordination](https://github.com/SuperInstance/categorical-coordination)** — Category theory for coordination

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project.
