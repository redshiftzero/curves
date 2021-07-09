use std::cmp::PartialEq;
use std::ops::{Add, Mul, Neg};

/// Weierstrass normal form - arithmetic in affine space over the reals

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Scalar {
    num: i32,
}

impl Scalar {
    fn new(num: i32) -> Self {
        Self { num }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct IdentityPoint {
    curve: WeierstrassNormalCurve,
}

impl IdentityPoint {
    fn new(curve: &WeierstrassNormalCurve) -> Self {
        Self { curve: *curve }
    }
}

/// P + O = P
impl<'a> Add<&'a AffinePoint> for IdentityPoint {
    type Output = &'a AffinePoint;

    fn add(self, other: &'a AffinePoint) -> &'a AffinePoint {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        other
    }
}

/// nP
impl Mul<Scalar> for IdentityPoint {
    type Output = IdentityPoint;

    fn mul(self, _scalar: Scalar) -> IdentityPoint {
        self
    }
}

/// nP
impl Mul<IdentityPoint> for Scalar {
    type Output = IdentityPoint;

    fn mul(self, point: IdentityPoint) -> IdentityPoint {
        point
    }
}

/// -O = O
impl Neg for IdentityPoint {
    type Output = IdentityPoint;

    fn neg(self) -> IdentityPoint {
        self
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct AffinePoint {
    x: i32,
    y: i32,
    curve: WeierstrassNormalCurve,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Point {
    AffinePoint(AffinePoint),     // x, y in affine space
    IdentityPoint(IdentityPoint), // point at infinity
}

impl AffinePoint {
    fn new(x: i32, y: i32, curve: &WeierstrassNormalCurve) -> Self {
        // Every time we create a point, we check if it's on the curve.
        // If not, we panic.
        let lhs = y.pow(2);
        let rhs = x.pow(3) + curve.a * x + curve.b;
        if rhs != lhs {
            panic!("err: Point not on curve! ({}, {})", x, y)
        }

        Self {
            x,
            y,
            curve: *curve,
        }
    }
}

/// O + P = P
impl Add<IdentityPoint> for AffinePoint {
    type Output = AffinePoint;

    fn add(self, other: IdentityPoint) -> AffinePoint {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        self
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;

    fn neg(self) -> AffinePoint {
        AffinePoint::new(self.x, -self.y, &self.curve)
    }
}

impl Add<AffinePoint> for AffinePoint {
    type Output = Point;

    fn add(self, other: AffinePoint) -> Point {
        if self.curve != other.curve {
            panic!("err: Adding points not on same curve")
        }

        let m: i32;
        if self.x == other.x && self.y == -other.y {
            return Point::IdentityPoint(IdentityPoint::new(&self.curve));
        } else if self.x == other.x && self.y == other.y {
            // P + P
            m = (3 * self.x.pow(2) + self.curve.a) / (2 * self.y);
        } else {
            // P + Q
            m = (self.y - other.y) / (self.x - other.x);
        }

        let r_x = m.pow(2) - self.x - other.x;
        let r_y = self.y + m * (r_x - self.x);

        // This will panic if the point ends up not being on the curve
        Point::AffinePoint(-AffinePoint::new(r_x, r_y, &self.curve))
    }
}

/// Curve in the form $y^2 = x^3 + ax + b$.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct WeierstrassNormalCurve {
    pub a: i32,
    pub b: i32,
}

impl WeierstrassNormalCurve {
    fn new(a: i32, b: i32) -> Self {
        if 4 * a.pow(3) + 27 * b.pow(2) == 0 {
            panic!("err: Curve is singular!")
        }

        Self { a, b }
    }
}

/// nP
impl Mul<Scalar> for AffinePoint {
    type Output = Point;

    fn mul(self, scalar: Scalar) -> Point {
        scalar * self
    }
}

/// nP
impl Mul<AffinePoint> for Scalar {
    type Output = Point;

    fn mul(self, point: AffinePoint) -> Point {
        if self.num == 0 {
            return Point::IdentityPoint(IdentityPoint::new(&point.curve));
        } else if self.num < 0 {
            self * -point
        } else {
            let mut q = point.clone();
            let mut r = point.clone();

            let mut i = 2;
            while i <= self.num {
                q = match q + q {
                    Point::IdentityPoint(identity_point) => {
                        return Point::IdentityPoint(identity_point)
                    }
                    Point::AffinePoint(affine_point) => affine_point,
                };

                if (self.num & i) == i {
                    r = match q + r {
                        Point::IdentityPoint(identity_point) => {
                            return Point::IdentityPoint(identity_point)
                        }
                        Point::AffinePoint(affine_point) => affine_point,
                    };
                }

                i = i << 1;
            }
            Point::AffinePoint(r)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "err: Point not on curve!")]
    fn cannot_create_point_not_on_curve() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let _point = AffinePoint::new(69, 420, &curve);
    }

    #[test]
    #[should_panic(expected = "err: Adding points not on same curve")]
    fn cannot_add_points_on_different_curves() {
        let curve1 = WeierstrassNormalCurve::new(-7, 10);
        let curve2 = WeierstrassNormalCurve::new(-7, 11);
        let point = AffinePoint::new(-3, 2, &curve1);
        let identity = IdentityPoint::new(&curve2);
        let _result = point + identity;
    }

    #[test]
    fn adding_point_to_identity_returns_point() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let point = AffinePoint::new(-3, 2, &curve);
        let identity = IdentityPoint::new(&curve);
        let result = point + identity;
        assert_eq!(result == point, true);
    }

    #[test]
    fn scalar_times_identity_is_identity() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let identity = IdentityPoint::new(&curve);
        let scalar = Scalar::new(4);
        let result = identity * scalar;
        assert_eq!(identity, result);
    }

    #[test]
    fn identity_times_scalar_is_identity() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let identity = IdentityPoint::new(&curve);
        let scalar = Scalar::new(4);
        let result = scalar * identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn negative_identity_is_identity() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let identity = IdentityPoint::new(&curve);
        let result = -identity;
        assert_eq!(identity, result);
    }

    #[test]
    fn negative_point_is_reflected_about_x() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let point = AffinePoint::new(-3, 2, &curve);
        let result = -point;
        assert_eq!(point.x, result.x);
        assert_eq!(point.y, -result.y);
    }

    // P + -P = O
    #[test]
    fn adding_point_to_neg_point_returns_identity() {
        let curve = WeierstrassNormalCurve::new(-7, 10);

        let p1 = AffinePoint::new(1, 2, &curve);
        let p2 = AffinePoint::new(1, -2, &curve);
        let expected = Point::IdentityPoint(IdentityPoint::new(&curve));
        let result = p1 + p2;
        assert_eq!(result, expected);
    }

    #[test]
    fn adding_point_to_point() {
        let curve = WeierstrassNormalCurve::new(-7, 10);

        let p1 = AffinePoint::new(1, 2, &curve);
        let p2 = AffinePoint::new(3, 4, &curve);
        let expected = Point::AffinePoint(AffinePoint::new(-3, 2, &curve));
        let result = p1 + p2;
        assert_eq!(result, expected);

        let p1 = AffinePoint::new(-1, 4, &curve);
        let p2 = AffinePoint::new(1, 2, &curve);
        let expected = Point::AffinePoint(AffinePoint::new(1, -2, &curve));
        let result = p1 + p2;
        assert_eq!(result, expected);
    }

    #[test]
    fn scalar_multiply() {
        let curve = WeierstrassNormalCurve::new(-7, 10);
        let p1 = AffinePoint::new(1, 2, &curve);
        let n = Scalar::new(3);
        let expected = Point::AffinePoint(AffinePoint::new(9, -26, &curve));
        let result = n * p1;
        assert_eq!(result, expected);

        let result = p1 * n;
        assert_eq!(result, expected);
    }
}
