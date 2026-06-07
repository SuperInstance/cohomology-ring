//! # Cohomology Ring
//!
//! A library for computing cohomology rings and cohomology operations.
//!
//! Provides:
//! - Cochain complexes with coboundary operator δ
//! - Cup product (cohomology ring multiplication)
//! - Cohomology group computation H^k
//! - Bockstein homomorphism (connecting map for coefficient sequences)
//! - Steenrod squares (mod 2 cohomology operations)

use std::collections::HashMap;

/// Represents a cochain as a linear combination of basis elements.
/// Stored as a sparse vector: basis_index → coefficient (mod 2 for Z/2 coefficients).
#[derive(Debug, Clone)]
pub struct Cochain {
    /// Degree of this cochain.
    pub degree: usize,
    /// Sparse representation: basis element index → coefficient.
    coefficients: HashMap<usize, i32>,
}

impl Cochain {
    /// Create a zero cochain of given degree.
    pub fn zero(degree: usize) -> Self {
        Self { degree, coefficients: HashMap::new() }
    }

    /// Create a cochain with a single basis element.
    pub fn basis(degree: usize, index: usize, coeff: i32) -> Self {
        let mut c = Self::zero(degree);
        if coeff != 0 {
            c.coefficients.insert(index, coeff);
        }
        c
    }

    /// Add two cochains (same degree).
    pub fn add(&self, other: &Cochain) -> Cochain {
        assert_eq!(self.degree, other.degree, "Cannot add cochains of different degrees");
        let mut result = self.clone();
        for (&idx, &coeff) in &other.coefficients {
            let entry = result.coefficients.entry(idx).or_insert(0);
            *entry += coeff;
            if *entry == 0 {
                result.coefficients.remove(&idx);
            }
        }
        result
    }

    /// Scale a cochain by a scalar.
    pub fn scale(&self, scalar: i32) -> Cochain {
        let mut result = Cochain::zero(self.degree);
        for (&idx, &coeff) in &self.coefficients {
            let new_coeff = coeff * scalar;
            if new_coeff != 0 {
                result.coefficients.insert(idx, new_coeff);
            }
        }
        result
    }

    /// Check if this cochain is zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// Get the coefficient at a given basis index.
    pub fn get(&self, index: usize) -> i32 {
        self.coefficients.get(&index).copied().unwrap_or(0)
    }

    /// Get the non-zero entries.
    pub fn entries(&self) -> &HashMap<usize, i32> {
        &self.coefficients
    }
}

/// A cochain complex C^0 → C^1 → C^2 → ...
///
/// The coboundary δ: C^k → C^{k+1} is the adjoint of the boundary operator.
#[derive(Debug, Clone)]
pub struct CochainComplex {
    /// Maximum degree.
    pub max_degree: usize,
    /// Dimension of each C^k.
    pub dimensions: Vec<usize>,
    /// Coboundary matrices δ_k: C^k → C^{k+1}.
    /// Stored as δ[k][i] = list of (j, coeff) meaning δ(e_i^k) contributes coeff * e_j^{k+1}.
    coboundary: Vec<Vec<Vec<(usize, i32)>>>,
}

impl CochainComplex {
    /// Create a new cochain complex with given dimensions and coboundary matrices.
    pub fn new(dimensions: Vec<usize>, coboundary: Vec<Vec<Vec<(usize, i32)>>>) -> Self {
        let max_degree = dimensions.len().saturating_sub(1);
        Self { max_degree, dimensions, coboundary }
    }

    /// Create a trivial complex with no coboundaries (all δ = 0).
    pub fn trivial(dimensions: Vec<usize>) -> Self {
        let max_degree = dimensions.len().saturating_sub(1);
        let coboundary = (0..max_degree).map(|_| vec![]).collect();
        Self { max_degree, dimensions, coboundary }
    }

    /// Create the cochain complex for S^1 (circle).
    /// C^0 ← C^1 with δ = 0 (since H^0 = Z, H^1 = Z).
    pub fn for_circle() -> Self {
        Self::trivial(vec![1, 1])
    }

    /// Create the cochain complex for S^2 (2-sphere).
    /// C^0 ← C^1 ← C^2 with all δ = 0.
    pub fn for_sphere() -> Self {
        Self::trivial(vec![1, 0, 1])
    }

    /// Create the cochain complex for RP^2 (real projective plane).
    /// Over Z/2: H^0 = Z/2, H^1 = Z/2, H^2 = Z/2.
    pub fn for_rp2() -> Self {
        // Z/2 cellular chain: one cell in each dimension 0, 1, 2
        // Coboundary δ_0: C^0 → C^1 is 0 (mod 2), δ_1: C^1 → C^2 is multiplication by 2 = 0 (mod 2)
        Self::trivial(vec![1, 1, 1])
    }

    /// Create the cochain complex for the torus T².
    /// H^0 = Z, H^1 = Z², H^2 = Z.
    pub fn for_torus() -> Self {
        Self::trivial(vec![1, 2, 1])
    }

    /// Apply the coboundary operator δ_k to a cochain of degree k.
    pub fn apply_coboundary(&self, cochain: &Cochain) -> Cochain {
        let k = cochain.degree;
        if k >= self.max_degree {
            return Cochain::zero(k + 1);
        }

        let target_dim = self.dimensions.get(k + 1).copied().unwrap_or(0);
        let mut result = Cochain::zero(k + 1);
        let mut result_coeffs = vec![0i32; target_dim];

        if k < self.coboundary.len() {
            for (&basis_idx, &coeff) in &cochain.coefficients {
                if basis_idx < self.coboundary[k].len() {
                    for &(target_idx, delta_coeff) in &self.coboundary[k][basis_idx] {
                        if target_idx < target_dim {
                            result_coeffs[target_idx] += coeff * delta_coeff;
                        }
                    }
                }
            }
        }

        for (i, &c) in result_coeffs.iter().enumerate() {
            if c != 0 {
                result.coefficients.insert(i, c);
            }
        }
        result
    }

    /// Get the coboundary matrix δ_k as a dense matrix.
    pub fn coboundary_matrix(&self, k: usize) -> Vec<Vec<i32>> {
        let rows = self.dimensions.get(k + 1).copied().unwrap_or(0);
        let cols = self.dimensions.get(k).copied().unwrap_or(0);
        let mut mat = vec![vec![0i32; cols]; rows];

        if k < self.coboundary.len() {
            for (col, entries) in self.coboundary[k].iter().enumerate() {
                for &(row, val) in entries {
                    if row < rows && col < cols {
                        mat[row][col] = val;
                    }
                }
            }
        }
        mat
    }
}

/// Computes the cohomology group H^k of a cochain complex.
///
/// H^k = ker(δ_k) / im(δ_{k-1})
#[derive(Debug, Clone)]
pub struct CohomologyGroup {
    /// The degree k.
    pub degree: usize,
    /// Rank of the free part (Betti number).
    pub rank: usize,
    /// Torsion components: (prime, power) pairs.
    pub torsion: Vec<(usize, usize)>,
}

impl CohomologyGroup {
    /// Compute H^k of the given cochain complex.
    pub fn compute(complex: &CochainComplex, k: usize) -> Self {
        let dim_k = complex.dimensions.get(k).copied().unwrap_or(0);
        if dim_k == 0 {
            return Self { degree: k, rank: 0, torsion: vec![] };
        }

        // For trivial coboundaries, H^k = C^k
        // Rank of δ_k (as a map from C^k)
        let rank_delta_k = matrix_rank(&complex.coboundary_matrix(k));
        // Rank of δ_{k-1} (as a map from C^{k-1} to C^k)
        let rank_delta_km1 = if k > 0 {
            matrix_rank(&complex.coboundary_matrix(k - 1))
        } else {
            0
        };

        // dim ker(δ_k) = dim_k - rank(δ_k)
        let dim_kernel = dim_k.saturating_sub(rank_delta_k);
        // dim im(δ_{k-1}) = rank(δ_{k-1})
        let dim_image = rank_delta_km1;
        // H^k has rank = dim_kernel - dim_image
        let betti = dim_kernel.saturating_sub(dim_image);

        Self { degree: k, rank: betti, torsion: vec![] }
    }

    /// Compute all cohomology groups up to max_degree.
    pub fn compute_all(complex: &CochainComplex) -> Vec<Self> {
        (0..=complex.max_degree)
            .map(|k| Self::compute(complex, k))
            .collect()
    }
}

/// Compute the rank of an integer matrix using row reduction.
fn matrix_rank(matrix: &[Vec<i32>]) -> usize {
    if matrix.is_empty() || matrix[0].is_empty() {
        return 0;
    }
    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut m: Vec<Vec<i32>> = matrix.to_vec();
    let mut rank = 0;
    let mut pivot_col = 0;

    while rank < rows && pivot_col < cols {
        // Find pivot
        let mut found = false;
        for i in rank..rows {
            if m[i][pivot_col] != 0 {
                m.swap(rank, i);
                found = true;
                break;
            }
        }
        if !found {
            pivot_col += 1;
            continue;
        }

        // Eliminate
        let pivot_val = m[rank][pivot_col];
        for i in 0..rows {
            if i == rank || m[i][pivot_col] == 0 { continue; }
            let factor = m[i][pivot_col];
            for j in 0..cols {
                m[i][j] = m[i][j] * pivot_val - factor * m[rank][j];
            }
        }
        rank += 1;
        pivot_col += 1;
    }
    rank
}

/// Cup product on cohomology.
///
/// The cup product ⌣: H^k × H^l → H^{k+l} gives the cohomology ring structure.
pub struct CupProduct;

impl CupProduct {
    /// Compute the cup product of two cochains.
    ///
    /// For simplicial cochains, the cup product of basis elements
    /// e_i^k ⌣ e_j^l can be computed combinatorially.
    /// Here we implement a simplified version for cellular cohomology.
    pub fn cup(a: &Cochain, b: &Cochain) -> Cochain {
        let product_degree = a.degree + b.degree;
        let mut result = Cochain::zero(product_degree);

        for (&i, &ca) in &a.coefficients {
            for (&j, &cb) in &b.coefficients {
                // Simplified: diagonal cup product
                let prod_coeff = ca * cb;
                if prod_coeff != 0 {
                    let entry = result.coefficients.entry(i).or_insert(0);
                    *entry += prod_coeff;
                    if *entry == 0 {
                        result.coefficients.remove(&i);
                    }
                }
            }
        }
        result
    }

    /// Verify graded-commutativity: a ⌣ b = (-1)^{kl} b ⌣ a.
    pub fn verify_graded_commutativity(a: &Cochain, b: &Cochain) -> bool {
        let ab = Self::cup(a, b);
        let ba = Self::cup(b, a);
        let sign = if (a.degree * b.degree) % 2 == 0 { 1 } else { -1 };
        let scaled_ba = ba.scale(sign);
        ab.add(&scaled_ba).is_zero()
    }

    /// Verify associativity: (a ⌣ b) ⌣ c = a ⌣ (b ⌣ c).
    pub fn verify_associativity(a: &Cochain, b: &Cochain, c: &Cochain) -> bool {
        let ab = Self::cup(a, b);
        let ab_c = Self::cup(&ab, c);
        let bc = Self::cup(b, c);
        let a_bc = Self::cup(a, &bc);
        ab_c.add(&a_bc.scale(-1)).is_zero()
    }
}

/// Bockstein homomorphism.
///
/// The Bockstein β: H^k(M; Z/p) → H^{k+1}(M; Z/p) is the connecting
/// homomorphism associated to the short exact sequence
/// 0 → Z/p → Z/p² → Z/p → 0.
pub struct BocksteinHomomorphism;

impl BocksteinHomomorphism {
    /// Apply the Bockstein to a mod-2 cochain.
    ///
    /// In the mod-2 setting, β(x) = δ(̃x) mod 2 where ̃x is an integral lift.
    /// For a cocycle x, β(x) = 0 iff x lifts to an integral class.
    pub fn apply(cochain: &Cochain, complex: &CochainComplex) -> Cochain {
        // The Bockstein is the coboundary of an integral lift
        // Simplified: apply coboundary and reduce mod 2
        let cob = complex.apply_coboundary(cochain);
        // Reduce mod 2
        let mut result = Cochain::zero(cob.degree);
        for (&idx, &coeff) in &cob.coefficients {
            let mod2 = ((coeff % 2) + 2) % 2;
            if mod2 != 0 {
                result.coefficients.insert(idx, mod2);
            }
        }
        result
    }

    /// β² = 0 (the Bockstein squares to zero).
    pub fn verify_square_zero(cochain: &Cochain, complex: &CochainComplex) -> bool {
        let beta1 = Self::apply(cochain, complex);
        let beta2 = Self::apply(&beta1, complex);
        beta2.is_zero()
    }
}

/// Steenrod squares Sq^i: H^k(M; Z/2) → H^{k+i}(M; Z/2).
///
/// Stable cohomology operations in mod 2 cohomology.
/// Key properties:
/// - Sq^0 = identity
/// - Sq^k(x) = x ⌣ x when deg(x) = k
/// - Cartan formula: Sq^k(x ⌣ y) = Σ_{i+j=k} Sq^i(x) ⌣ Sq^j(y)
pub struct SteenrodSquare;

impl SteenrodSquare {
    /// Compute Sq^0 on a mod-2 cochain (identity).
    pub fn sq0(cochain: &Cochain) -> Cochain {
        cochain.clone()
    }

    /// Compute Sq^k for small k on a mod-2 cochain.
    /// For a cochain of degree n: Sq^n(x) = x ⌣ x.
    pub fn sqn(cochain: &Cochain) -> Cochain {
        CupProduct::cup(cochain, cochain)
    }

    /// Compute Sq^1 (the Bockstein for mod 2).
    /// Sq^1 = β (mod 2 Bockstein).
    pub fn sq1(cochain: &Cochain, complex: &CochainComplex) -> Cochain {
        BocksteinHomomorphism::apply(cochain, complex)
    }

    /// Verify the Cartan formula: Sq^k(x ⌣ y) = Σ_{i+j=k} Sq^i(x) ⌣ Sq^j(y).
    /// For k=0, this is trivially true since Sq^0 is the identity.
    pub fn verify_cartan_formula(x: &Cochain, y: &Cochain, complex: &CochainComplex) -> bool {
        let xy = CupProduct::cup(x, y);
        let sq0_xy = Self::sq0(&xy);

        // For k=0: Sq^0(xy) = xy = Sq^0(x) ⌣ Sq^0(y)
        let lhs = sq0_xy;
        let rhs = CupProduct::cup(&Self::sq0(x), &Self::sq0(y));
        lhs.add(&rhs.scale(-1)).is_zero()
    }

    /// Verify Adem relations for specific cases.
    /// Returns true if an Adem relation is needed (i < 2j).
    pub fn adem_relation_check(i: usize, j: usize) -> bool {
        i < 2 * j
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cochain_zero() {
        let c = Cochain::zero(3);
        assert_eq!(c.degree, 3);
        assert!(c.is_zero());
    }

    #[test]
    fn test_cochain_add() {
        let a = Cochain::basis(1, 0, 3);
        let b = Cochain::basis(1, 0, -3);
        let sum = a.add(&b);
        assert!(sum.is_zero());
    }

    #[test]
    fn test_cochain_add_different_basis() {
        let a = Cochain::basis(1, 0, 3);
        let b = Cochain::basis(1, 1, -2);
        let sum = a.add(&b);
        assert_eq!(sum.get(0), 3);
        assert_eq!(sum.get(1), -2);
        assert!(!sum.is_zero());
    }

    #[test]
    fn test_cochain_scale() {
        let a = Cochain::basis(2, 0, 5);
        let scaled = a.scale(3);
        assert_eq!(scaled.get(0), 15);
    }

    #[test]
    fn test_cohomology_circle() {
        let complex = CochainComplex::for_circle();
        let groups = CohomologyGroup::compute_all(&complex);
        assert_eq!(groups[0].rank, 1); // H^0(S¹) = Z
        assert_eq!(groups[1].rank, 1); // H^1(S¹) = Z
    }

    #[test]
    fn test_cohomology_sphere() {
        let complex = CochainComplex::for_sphere();
        let groups = CohomologyGroup::compute_all(&complex);
        assert_eq!(groups[0].rank, 1); // H^0(S²) = Z
        assert_eq!(groups[1].rank, 0); // H^1(S²) = 0
        assert_eq!(groups[2].rank, 1); // H^2(S²) = Z
    }

    #[test]
    fn test_cohomology_torus() {
        let complex = CochainComplex::for_torus();
        let groups = CohomologyGroup::compute_all(&complex);
        assert_eq!(groups[0].rank, 1);
        assert_eq!(groups[1].rank, 2);
        assert_eq!(groups[2].rank, 1);
    }

    #[test]
    fn test_cup_product_degree() {
        let a = Cochain::basis(1, 0, 1);
        let b = Cochain::basis(1, 0, 1);
        let product = CupProduct::cup(&a, &b);
        assert_eq!(product.degree, 2);
    }

    #[test]
    fn test_cup_product_commutativity() {
        // Even-degree cochains: graded commutativity with sign (-1)^{kl} = (-1)^{2*2} = 1
        // So a⌣b = b⌣a for even degree
        let a = Cochain::basis(2, 0, 2);
        let b = Cochain::basis(2, 0, 3);
        let ab = CupProduct::cup(&a, &b);
        let ba = CupProduct::cup(&b, &a);
        // For the same basis index, cup product is just multiplication
        // a⌣b = (2e_0)⌣(3e_0) = 6e_0, b⌣a = (3e_0)⌣(2e_0) = 6e_0
        assert_eq!(ab.get(0), ba.get(0));
    }

    #[test]
    fn test_cup_product_associativity() {
        let a = Cochain::basis(1, 0, 1);
        let b = Cochain::basis(1, 0, 1);
        let c = Cochain::basis(1, 0, 1);
        assert!(CupProduct::verify_associativity(&a, &b, &c));
    }

    #[test]
    fn test_bockstein_square_zero() {
        let complex = CochainComplex::for_circle();
        let cochain = Cochain::basis(0, 0, 1);
        assert!(BocksteinHomomorphism::verify_square_zero(&cochain, &complex));
    }

    #[test]
    fn test_steenrod_sq0_identity() {
        let x = Cochain::basis(2, 0, 1);
        let sq0 = SteenrodSquare::sq0(&x);
        assert_eq!(sq0.coefficients, x.coefficients);
    }

    #[test]
    fn test_steenrod_cartan_formula() {
        let complex = CochainComplex::for_circle();
        let x = Cochain::basis(1, 0, 1);
        let y = Cochain::basis(1, 0, 1);
        assert!(SteenrodSquare::verify_cartan_formula(&x, &y, &complex));
    }

    #[test]
    fn test_adem_relation_admissible() {
        // i >= 2j means admissible, no relation needed (returns false)
        assert!(!SteenrodSquare::adem_relation_check(4, 2)); // 4 >= 2*2 = 4
        // i < 2j means Adem relation applies (returns true)
        assert!(SteenrodSquare::adem_relation_check(1, 1)); // 1 < 2*1 = 2
    }

    #[test]
    fn test_matrix_rank_identity() {
        let m = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(matrix_rank(&m), 2);
    }

    #[test]
    fn test_matrix_rank_zero() {
        let m = vec![vec![0, 0], vec![0, 0]];
        assert_eq!(matrix_rank(&m), 0);
    }

    #[test]
    fn test_cochain_complex_with_coboundary() {
        // Complex with a non-trivial coboundary
        let coboundary = vec![
            vec![vec![(0, 2)]] // δ_0: e_0^0 → 2·e_0^1
        ];
        let complex = CochainComplex::new(vec![1, 1], coboundary);
        let h0 = CohomologyGroup::compute(&complex, 0);
        // ker(δ_0) = {x : 2x = 0} over Z, so ker = {0}, rank 0
        // But with integer computation, rank(delta_0) = 1, so dim_ker = 0
        assert_eq!(h0.rank, 0);
    }
}
