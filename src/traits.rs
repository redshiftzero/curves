pub trait EllipticCurve<T> {
    fn is_point_on_curve(&self, p: &T) -> bool;
}
