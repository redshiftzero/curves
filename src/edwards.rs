use std::ops::{Add, Mul, Neg};

use crate::field::PrimeFieldElement19;
use crate::scalar::Scalar;
use crate::traits::EllipticCurve;

/// Curve in the form $x**2 + y**2 = c**2 ( 1 + d * x**2 * y**2 )$.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct EdwardsCurve {
    pub c: PrimeFieldElement19,
    pub d: PrimeFieldElement19,
}

/// A point on a curve in Edwards form over $GF(19)$ in affine space.
///
/// When instantiating `EdwardsAffinePoint` via `new()` we check that
/// the point is on the curve.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct EdwardsAffinePoint {
    x: PrimeFieldElement19,
    y: PrimeFieldElement19,
    curve: EdwardsCurve,
}

impl EdwardsCurve {
    fn new(c: PrimeFieldElement19, d: PrimeFieldElement19) -> Self {
        Self { c, d }
    }
}

impl EllipticCurve<EdwardsAffinePoint> for EdwardsCurve {
    fn is_point_on_curve(&self, p: &EdwardsAffinePoint) -> bool {
        let lhs = p.x.pow(2) + p.y.pow(2);
        let rhs = self.c.pow(2) * (PrimeFieldElement19::new(1) + self.d * p.x.pow(2) * p.y.pow(2));
        rhs == lhs
    }

    /// Edwards neutral element at $(0,c)$.
    fn identity(&self) -> EdwardsAffinePoint {
        EdwardsAffinePoint::new(PrimeFieldElement19::new(0), self.c, &self)
    }
}

impl EdwardsAffinePoint {
    fn new(x: PrimeFieldElement19, y: PrimeFieldElement19, curve: &EdwardsCurve) -> Self {
        let point = Self {
            x,
            y,
            curve: *curve,
        };

        // Every time we create a point, we check if it's on the curve.
        // If not, we panic.
        if !curve.is_point_on_curve(&point) {
            panic!("err: Point not on curve! ({:?}, {:?})", x, y)
        }

        point
    }
}

/// Edwards negative of (x, y) is (-x, y).
///
/// From the elliptic curve equation we have:
/// $ x^2 + y^2 = c^2 + d c^2 x^2 y^2 $
///
/// Rearranging to get x^2 on one side:
/// $ x^2 - d c^2 x^2 y^2 = c^2 - y^2 $
///
/// Factoring out x^2 from the LHS:
/// $ x^2 ( 1 - d c^2 y^2 ) = c^2 - y^2 $
///
/// Giving us finally:
/// $ x^2 = \frac{c^2 - y^2}{1 - d c^2 y^2} $
///
/// Taking the square root will yield us the negation of x.
impl Neg for EdwardsAffinePoint {
    type Output = EdwardsAffinePoint;

    fn neg(self) -> EdwardsAffinePoint {
        // Recalling that our field inverse is implemented using the Neg trait:
        let x_squared = (self.curve.c.pow(2) - self.y.pow(2))
            * -(PrimeFieldElement19::new(1) - self.curve.d * self.curve.c.pow(2) * self.y.pow(2));

        // The square root of the RHS is the new y
        // Special case for GF(19) since 19 is congruent to 3 mod 4
        let m = 19;
        let first_root = x_squared.pow((m + 1) / 4);
        let second_root = PrimeFieldElement19::new(19 - first_root.x);

        if first_root == self.x {
            EdwardsAffinePoint::new(second_root, self.y, &self.curve)
        } else {
            EdwardsAffinePoint::new(first_root, self.y, &self.curve)
        }
    }
}

impl Add<EdwardsAffinePoint> for EdwardsAffinePoint {
    type Output = EdwardsAffinePoint;

    fn add(self, other: EdwardsAffinePoint) -> EdwardsAffinePoint {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        let r_x = (self.x * other.y + self.y * other.x)
            * -(self.curve.c
                * (PrimeFieldElement19::new(1)
                    + self.curve.d * self.x * other.x * self.y * other.y));
        let r_y = (self.y * other.y - self.x * other.x)
            * -(self.curve.c
                * (PrimeFieldElement19::new(1)
                    - self.curve.d * self.x * other.x * self.y * other.y));

        // This will panic if the point ends up not being on the curve
        EdwardsAffinePoint::new(r_x, r_y, &self.curve)
    }
}

impl Mul<Scalar> for EdwardsAffinePoint {
    type Output = EdwardsAffinePoint;

    fn mul(self, scalar: Scalar) -> EdwardsAffinePoint {
        scalar * self
    }
}

impl Mul<EdwardsAffinePoint> for Scalar {
    type Output = EdwardsAffinePoint;

    fn mul(self, point: EdwardsAffinePoint) -> EdwardsAffinePoint {
        if self.num == 0 {
            return point.curve.identity();
        } else if self.num < 0 {
            -self * -point
        } else {
            let mut q = point;
            let mut r: EdwardsAffinePoint;
            if self.num & 1 == 1 {
                r = point;
            } else {
                r = point.curve.identity();
            }

            let mut i = 2;

            while i <= self.num {
                q = q + q;

                if (self.num & i) == i {
                    r = q + r;
                }

                i = i << 1;
            }
            r
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "err: Point not on curve!")]
    fn cannot_create_point_not_on_curve() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let _point = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(69),
            PrimeFieldElement19::new(420),
            &curve,
        );
    }

    #[test]
    fn identity_point_on_curve() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let _point = curve.identity();
    }

    #[test]
    fn negative_identity_is_identity() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let result = -identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn point_negation() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let p = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(2),
            PrimeFieldElement19::new(0),
            &curve,
        );
        // Will panic if the result is not on the curve!
        let result = -p;
        assert_ne!(result, p);
        assert_eq!(result.x, PrimeFieldElement19::new(17)); // -2 mod 19
        assert_eq!(result.y, p.y);
    }

    #[test]
    #[should_panic(expected = "err: Adding points not on same curve")]
    fn cannot_add_points_on_different_curves() {
        let curve1 = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let point1 = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(2),
            PrimeFieldElement19::new(0),
            &curve1,
        );
        let curve2 = EdwardsCurve::new(PrimeFieldElement19::new(3), PrimeFieldElement19::new(3));
        let point2 = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(3),
            PrimeFieldElement19::new(0),
            &curve2,
        );
        let _result = point1 + point2;
    }

    #[test]
    fn adding_point_to_identity_returns_point() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let point = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(2),
            PrimeFieldElement19::new(0),
            &curve,
        );
        let identity = curve.identity();
        let result = point + identity;
        assert_eq!(result, point);

        let result = identity + point;
        assert_eq!(result, point);
    }

    #[test]
    fn adding_point_to_neg_point_returns_identity() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let p1 = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(2),
            PrimeFieldElement19::new(0),
            &curve,
        );
        let p2 = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(17),
            PrimeFieldElement19::new(0),
            &curve,
        );
        let identity = curve.identity();
        let result = p1 + p2;
        assert_eq!(result, identity);
    }

    #[test]
    fn point_doubling() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let p = EdwardsAffinePoint::new(
            PrimeFieldElement19::new(2),
            PrimeFieldElement19::new(0),
            &curve,
        );
        let expected = p + p;

        let scalar = Scalar::new(2);
        let result = scalar * p;

        assert_eq!(result, expected);
    }

    #[test]
    fn scalar_times_identity_is_identity() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let scalar = Scalar::new(4);
        let result = identity * scalar;
        assert_eq!(identity, result);
    }

    #[test]
    fn identity_times_scalar_is_identity() {
        let curve = EdwardsCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(2));
        let identity = curve.identity();
        let scalar = Scalar::new(4);
        let result = scalar * identity;
        assert_eq!(identity, result);
    }
}
