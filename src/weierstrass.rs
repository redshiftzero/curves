use std::ops::{Add, Mul, Neg};

use crate::field::PrimeFieldElement19;
use crate::scalar::Scalar;
use crate::traits::EllipticCurve;

/// A point on a curve in Weierstrass normal form over $GF(19)$ in affine space.
///
/// Note: If Option is None, then this point is at infinity.
///
/// When instantiating `WeierstrassAffinePoint` via `new()`, we check that:
/// 1. The point is either at infinity, meaning x and y are _both_ `None`, and
/// 2. The point is actually on the curve.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct WeierstrassAffinePoint {
    x: Option<PrimeFieldElement19>,
    y: Option<PrimeFieldElement19>,
    curve: WeierstrassNormalCurve,
}

/// Curve in the form $y^2 = x^3 + ax + b$.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct WeierstrassNormalCurve {
    pub a: PrimeFieldElement19,
    pub b: PrimeFieldElement19,
}

impl WeierstrassNormalCurve {
    fn new(a: PrimeFieldElement19, b: PrimeFieldElement19) -> Self {
        if (4 * a.pow(3) + 27 * b.pow(2)).x == 0 {
            panic!("err: Curve is singular!")
        }

        Self { a, b }
    }
}

impl EllipticCurve<WeierstrassAffinePoint> for WeierstrassNormalCurve {
    fn is_point_on_curve(&self, p: &WeierstrassAffinePoint) -> bool {
        match p.x {
            None => return true, // At infinity
            Some(x) => {
                // Safe unwrap because we ensure that points either have Some
                // for both fields x, y or both Nones in new().
                let y = p.y.unwrap();

                let lhs = y.pow(2);
                let rhs = x.pow(3) + self.a * x + self.b;

                if rhs != lhs {
                    return false;
                }
            }
        };

        true
    }
}

impl Neg for WeierstrassAffinePoint {
    type Output = WeierstrassAffinePoint;

    fn neg(self) -> WeierstrassAffinePoint {
        match self.x {
            None => return self, // -O = O
            Some(x) => {
                let y = self.y.unwrap();

                // Given an (x,y) pair we want to know: What is the other point on the curve
                // that has the same x-coordinate?
                let rhs = x.pow(3) + self.curve.a * x + self.curve.b;

                // The square root of the RHS is the new y
                // Special case for GF(19) since 19 is congruent to 3 mod 4
                let m = 19;
                let first_root = rhs.pow((m + 1) / 4);
                let second_root = PrimeFieldElement19::new(19 - first_root.x);

                if first_root == y {
                    WeierstrassAffinePoint::new(Some(x), Some(second_root), &self.curve)
                } else {
                    WeierstrassAffinePoint::new(Some(x), Some(first_root), &self.curve)
                }
            }
        }
    }
}

impl Add<WeierstrassAffinePoint> for WeierstrassAffinePoint {
    type Output = WeierstrassAffinePoint;

    fn add(self, other: WeierstrassAffinePoint) -> WeierstrassAffinePoint {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        if self.x.is_none() {
            return other; // O + P = P
        } else if other.x.is_none() {
            return self; // P + O = P
        } else {
            // Either P + P, P + Q, or P + -P
            let m: PrimeFieldElement19;
            let r_x: PrimeFieldElement19;

            let x = self.x.unwrap();
            let y = self.y.unwrap();
            let other_x = other.x.unwrap();
            let other_y = other.y.unwrap();

            if self == -other {
                // -P + P
                return WeierstrassAffinePoint::new(None, None, &self.curve);
            } else if x == other_x && y == other_y {
                // P + P
                // The below line is: m = (3 * x.pow(2) + self.curve.a) / (2 * y);
                // using the fact that mod inverse is implemented using the Neg trait
                m = (3 * x.pow(2) + self.curve.a) * -(2 * y);
                r_x = m.pow(2) - 2 * x;
            } else {
                // P + Q
                // The below line is: m = (y - other_y) / (x - other_x)
                // using the fact that mod inverse is implemented using the Neg trait
                m = (y - other_y) * -(x - other_x);
                r_x = m.pow(2) - x - other_x;
            }

            let r_y = y + m * (r_x - x);

            // This will panic if the point ends up not being on the curve
            -WeierstrassAffinePoint::new(Some(r_x), Some(r_y), &self.curve)
        }
    }
}

impl WeierstrassAffinePoint {
    fn new(
        x: Option<PrimeFieldElement19>,
        y: Option<PrimeFieldElement19>,
        curve: &WeierstrassNormalCurve,
    ) -> Self {
        let point = Self {
            x,
            y,
            curve: *curve,
        };

        if (x.is_none() & !y.is_none()) | (!x.is_none() & y.is_none()) {
            panic!("err: both x, y must be None for a valid point at infinity")
        }

        // Every time we create a point, we check if it's on the curve.
        // If not, we panic.
        if !curve.is_point_on_curve(&point) {
            panic!("err: Point not on curve! ({:?}, {:?})", x, y)
        }

        point
    }
}

impl Mul<Scalar> for WeierstrassAffinePoint {
    type Output = WeierstrassAffinePoint;

    fn mul(self, scalar: Scalar) -> WeierstrassAffinePoint {
        scalar * self
    }
}

impl Mul<WeierstrassAffinePoint> for Scalar {
    type Output = WeierstrassAffinePoint;

    fn mul(self, point: WeierstrassAffinePoint) -> WeierstrassAffinePoint {
        if self.num == 0 {
            return WeierstrassAffinePoint::new(None, None, &point.curve);
        } else if self.num < 0 {
            -self * -point
        } else {
            let mut q = point;
            let mut r: WeierstrassAffinePoint;
            if self.num & 1 == 1 {
                r = point;
            } else {
                r = WeierstrassAffinePoint::new(None, None, &point.curve);
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
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(7), PrimeFieldElement19::new(10));
        let _point = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(69)),
            Some(PrimeFieldElement19::new(420)),
            &curve,
        );
    }

    #[test]
    #[should_panic(expected = "err: Adding points not on same curve")]
    fn cannot_add_points_on_different_curves() {
        let curve1 =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(7), PrimeFieldElement19::new(10));
        let curve2 =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(7), PrimeFieldElement19::new(11));
        let point1 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(12)),
            Some(PrimeFieldElement19::new(13)),
            &curve1,
        );
        let point2 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(10)),
            Some(PrimeFieldElement19::new(13)),
            &curve2,
        );
        let _result = point1 + point2;
    }

    #[test]
    fn adding_point_to_identity_returns_point() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(7), PrimeFieldElement19::new(10));
        let point = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(12)),
            Some(PrimeFieldElement19::new(13)),
            &curve,
        );
        let identity = WeierstrassAffinePoint::new(None, None, &curve);
        let result = point + identity;
        assert_eq!(result, point);
    }

    #[test]
    fn add_point_to_point() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p1 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(15)),
            Some(PrimeFieldElement19::new(11)),
            &curve,
        );
        let p2 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(14)),
            Some(PrimeFieldElement19::new(1)),
            &curve,
        );
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(14)),
            Some(PrimeFieldElement19::new(18)),
            &curve,
        );
        let result = p1 + p2;
        assert_eq!(result, expected);
    }

    #[test]
    fn point_doubling() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p1 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(15)),
            Some(PrimeFieldElement19::new(11)),
            &curve,
        );
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(5)),
            Some(PrimeFieldElement19::new(10)),
            &curve,
        );
        let result = p1 + p1;
        assert_eq!(result, expected);
    }

    #[test]
    fn point_negation() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(15)),
            Some(PrimeFieldElement19::new(11)),
            &curve,
        );
        // Will panic if the result is not on the curve!
        let result = -p;
        assert_ne!(result, p);
    }

    #[test]
    fn negative_identity_is_identity() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let identity = WeierstrassAffinePoint::new(None, None, &curve);
        let result = -identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn adding_point_to_neg_point_returns_identity() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p1 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(5)),
            Some(PrimeFieldElement19::new(9)),
            &curve,
        );
        let p2 = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(5)),
            Some(PrimeFieldElement19::new(10)),
            &curve,
        );
        let identity = WeierstrassAffinePoint::new(None, None, &curve);
        let result = p1 + p2;
        assert_eq!(result, identity);
    }

    #[test]
    fn scalar_times_identity_is_identity() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let identity = WeierstrassAffinePoint::new(None, None, &curve);
        let scalar = Scalar::new(4);
        let result = identity * scalar;
        assert_eq!(identity, result);
    }

    #[test]
    fn identity_times_scalar_is_identity() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let identity = WeierstrassAffinePoint::new(None, None, &curve);
        let scalar = Scalar::new(4);
        let result = scalar * identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn scalar_multiply_doubling() {
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(15)),
            Some(PrimeFieldElement19::new(11)),
            &curve,
        );
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(5)),
            Some(PrimeFieldElement19::new(10)),
            &curve,
        );
        let scalar = Scalar::new(2);
        let result = scalar * p;
        assert_eq!(expected, result);
    }

    #[test]
    fn scalar_multiply() {
        // 3P
        let curve =
            WeierstrassNormalCurve::new(PrimeFieldElement19::new(2), PrimeFieldElement19::new(3));
        let p = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(3)),
            Some(PrimeFieldElement19::new(6)),
            &curve,
        );
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(15)),
            Some(PrimeFieldElement19::new(11)),
            &curve,
        );
        let scalar = Scalar::new(3);
        let result = scalar * p;
        assert_eq!(expected, result);

        // 4P
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(5)),
            Some(PrimeFieldElement19::new(9)),
            &curve,
        );
        let scalar = Scalar::new(4);
        let result = scalar * p;
        assert_eq!(expected, result);

        // 5P
        let expected = WeierstrassAffinePoint::new(
            Some(PrimeFieldElement19::new(18)),
            Some(PrimeFieldElement19::new(0)),
            &curve,
        );
        let scalar = Scalar::new(5);
        let result = scalar * p;
        assert_eq!(expected, result);
    }
}
