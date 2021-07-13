use std::ops::Neg;

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
