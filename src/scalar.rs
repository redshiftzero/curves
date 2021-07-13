#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Scalar {
    pub num: i32,
}

impl Scalar {
    pub fn new(num: i32) -> Self {
        Self { num }
    }
}
