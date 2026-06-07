# cohomology-ring

> **Cup product. Cohomology operations. The ring structure of H*(X).**

[![crates.io](https://img.shields.io/crates/v/cohomology-ring.svg)](https://crates.io/crates/cohomology-ring)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Computes cohomology groups, cup products, and cohomology operations including the Bockstein homomorphism and basic Steenrod squares. The cup product gives cohomology its ring structure — turning it from a graded vector space into a powerful algebraic invariant.

## The Cup Product

Cohomology isn't just groups H⁰, H¹, H², ... — they form a **graded ring** via the cup product:

∪: Hᵏ × Hˡ → Hᵏ⁺ˡ

This ring structure distinguishes spaces that homology alone cannot. For example, CP² and S²∨S⁴ have the same homology groups but different cohomology rings.

## Operations

- **Cup product**: the multiplicative structure on cohomology
- **Bockstein**: connecting homomorphism β: Hᵏ(X; Z/p) → Hᵏ⁺¹(X; Z/p)
- **Steenrod squares**: Sqⁱ: Hᵏ(X; Z/2) → Hᵏ⁺ⁱ(X; Z/2)

## License

MIT © [SuperInstance](https://github.com/SuperInstance)
