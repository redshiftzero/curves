use std::ops::{Mul, Neg};

use crate::traits::AffinePoint;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Scalar {
    pub num: i32,
}

impl Scalar {
    pub fn new(num: i32) -> Self {
        Self { num }
    }
}

impl Neg for Scalar {
    type Output = Scalar;

    fn neg(self) -> Scalar {
        Scalar::new(-1 * self.num)
    }
}

/// Default implementation of scalar multiplication using the double-and-add method
/// for affine points that are Copy and implement both the Add and Neg traits.
impl<T: AffinePoint + Copy + std::ops::Add<Output = T> + std::ops::Neg<Output = T>> Mul<T>
    for Scalar
{
    type Output = T;

    fn mul(self, point: T) -> T {
        if self.num == 0 {
            return point.identity();
        } else if self.num < 0 {
            -self * -point
        } else {
            let mut q = point;
            let mut r: T;
            if self.num & 1 == 1 {
                r = point;
            } else {
                r = point.identity();
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
