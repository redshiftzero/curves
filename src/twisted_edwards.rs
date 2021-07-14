use std::ops::{Add, Mul, Neg};

use crate::field::PrimeFieldElement19;
use crate::scalar::Scalar;
use crate::traits::{AffinePoint, EllipticCurve};

/// Curve in the form $a * x**2 + y**2 = 1 + d * x**2 * y**2$.
///
/// The coefficient $a$ is the "twist". When $a=1$ the curve is
/// a regular `EdwardsCurve`.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct TwistedEdwardsCurve {
    pub a: PrimeFieldElement19,
    pub d: PrimeFieldElement19,
}

/// A point on a curve in Twisted Edwards form over $GF(19)$ in affine space.
///
/// When instantiating `TwistedEdwardsAffinePoint` via `new()` we check that
/// the point is on the curve.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct TwistedEdwardsAffinePoint {
    x: PrimeFieldElement19,
    y: PrimeFieldElement19,
    curve: TwistedEdwardsCurve,
}

impl AffinePoint for TwistedEdwardsAffinePoint {
    fn identity(&self) -> TwistedEdwardsAffinePoint {
        self.curve.identity()
    }
}

impl TwistedEdwardsCurve {
    fn new(a: PrimeFieldElement19, d: PrimeFieldElement19) -> Self {
        Self { a, d }
    }
}

impl EllipticCurve<TwistedEdwardsAffinePoint> for TwistedEdwardsCurve {
    fn is_point_on_curve(&self, p: &TwistedEdwardsAffinePoint) -> bool {
        let lhs = self.a * p.x.pow(2) + p.y.pow(2);
        let rhs = PrimeFieldElement19::new(1) + self.d * p.x.pow(2) * p.y.pow(2);
        rhs == lhs
    }

    /// Twisted Edwards neutral element is at $(0,1)$.
    fn identity(&self) -> TwistedEdwardsAffinePoint {
        TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(1),
            &self,
        )
    }
}

impl TwistedEdwardsAffinePoint {
    fn new(x: PrimeFieldElement19, y: PrimeFieldElement19, curve: &TwistedEdwardsCurve) -> Self {
        let point = Self {
            x,
            y,
            curve: *curve,
        };

        if !curve.is_point_on_curve(&point) {
            panic!("err: Point not on curve! ({:?}, {:?})", x, y)
        }

        point
    }
}

/// Twisted Edwards negative of (x, y) is (-x, y).
///
/// From the elliptic curve equation we have:
/// $ a x^2 + y^2 = 1 + d x^2 y^2 $
///
/// Rearranging to get x^2 on one side:
/// $ a x^2 - d x^2 y^2 = 1 - y^2 $
///
/// Factoring out x^2 from the LHS:
/// $ x^2 ( a - d y^2 ) = 1 - y^2 $
///
/// Giving us finally:
/// $ x^2 = \frac{1 - y^2}{a - d y^2} $
///
/// Taking the square root will yield us the negation of x.
impl Neg for TwistedEdwardsAffinePoint {
    type Output = TwistedEdwardsAffinePoint;

    fn neg(self) -> TwistedEdwardsAffinePoint {
        if self == self.identity() {
            return self;
        }

        // Recalling that our field inverse is implemented using the Neg trait:
        let x_squared = (PrimeFieldElement19::new(1) - self.y.pow(2))
            * -(self.curve.a - self.curve.d * self.y.pow(2));

        // The square root of the RHS is the new y
        // Special case for GF(19) since 19 is congruent to 3 mod 4
        let m = 19;
        let first_root = x_squared.pow((m + 1) / 4);
        let second_root = PrimeFieldElement19::new(19 - first_root.x);

        if first_root == self.x {
            TwistedEdwardsAffinePoint::new(second_root, self.y, &self.curve)
        } else {
            TwistedEdwardsAffinePoint::new(first_root, self.y, &self.curve)
        }
    }
}

impl Add<TwistedEdwardsAffinePoint> for TwistedEdwardsAffinePoint {
    type Output = TwistedEdwardsAffinePoint;

    fn add(self, other: TwistedEdwardsAffinePoint) -> TwistedEdwardsAffinePoint {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        let r_x = (self.x * other.y + self.y * other.x)
            * -(PrimeFieldElement19::new(1) + self.curve.d * self.x * other.x * self.y * other.y);
        let r_y = (self.y * other.y - self.curve.a * self.x * other.x)
            * (PrimeFieldElement19::new(1) - self.curve.d * self.x * other.x * self.y * other.y);

        // This will panic if the point ends up not being on the curve
        TwistedEdwardsAffinePoint::new(r_x, r_y, &self.curve)
    }
}

impl Mul<Scalar> for TwistedEdwardsAffinePoint {
    type Output = TwistedEdwardsAffinePoint;

    fn mul(self, scalar: Scalar) -> TwistedEdwardsAffinePoint {
        scalar * self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "err: Point not on curve!")]
    fn cannot_create_point_not_on_curve() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let _point = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(69),
            PrimeFieldElement19::new(420),
            &curve,
        );
    }

    #[test]
    fn identity_point_on_curve() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let _point = curve.identity();
    }

    #[test]
    fn negative_identity_is_identity() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let result = -identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn point_negation() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(8), PrimeFieldElement19::new(4));
        let p = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(18),
            &curve,
        );
        let result = -p;
        assert_eq!(result.y, p.y);
    }

    #[test]
    #[should_panic(expected = "err: Adding points not on same curve")]
    fn cannot_add_points_on_different_curves() {
        let curve1 =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let point1 = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(18),
            &curve1,
        );
        let curve2 =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(3), PrimeFieldElement19::new(3));
        let point2 = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(1),
            &curve2,
        );
        let _result = point1 + point2;
    }

    #[test]
    fn adding_point_to_identity_returns_point() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let point = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(18),
            &curve,
        );
        let identity = curve.identity();
        let result = point + identity;
        assert_eq!(result, point);

        let result = identity + point;
        assert_eq!(result, point);
    }

    #[test]
    fn point_doubling() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let p = TwistedEdwardsAffinePoint::new(
            PrimeFieldElement19::new(0),
            PrimeFieldElement19::new(18),
            &curve,
        );
        let expected = p + p;

        let scalar = Scalar::new(2);
        let result = scalar * p;

        assert_eq!(result, expected);
    }

    #[test]
    fn scalar_times_identity_is_identity() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let scalar = Scalar::new(4);
        let result = identity * scalar;
        assert_eq!(identity, result);
    }

    #[test]
    fn identity_times_scalar_is_identity() {
        let curve =
            TwistedEdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let scalar = Scalar::new(4);
        let result = scalar * identity;
        assert_eq!(identity, result);
    }
}
