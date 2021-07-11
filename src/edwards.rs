use crate::traits::EllipticCurve;

// Curve in the form $x**2 + y**2 = 1 + d * x**2 * y**2$.
pub struct EdwardsCurve {
    pub d: i32,  //TODO be field element
}

