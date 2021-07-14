/// Trait defining an elliptic curve over a generic point type T.
pub trait EllipticCurve<T> {
    fn is_point_on_curve(&self, p: &T) -> bool;
    fn identity(&self) -> T;
}

/// Trait defining an affine point.
pub trait AffinePoint {
    /// The `identity` method should return the identity point on the
    /// associated elliptic curve.
    ///
    /// This is a bit of an odd place for this but it's required for the
    /// default scalar multiplication implementation (see src/scalar.rs).
    fn identity(&self) -> Self;
}
