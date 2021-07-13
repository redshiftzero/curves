use core::ops::{Add, Mul, Neg, Sub};

const MODULUS: u64 = 19;

/// `PrimeFieldElement19` represents an element of $GF(19)$ or $F_{19}$
/// i.e. a finite field of prime characteristic 19.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PrimeFieldElement19 {
    pub x: u64,
}

impl PrimeFieldElement19 {
    pub fn new(x: u64) -> Self {
        PrimeFieldElement19 {
            x: x.rem_euclid(MODULUS),
        }
    }

    pub fn pow(&self, pow: u32) -> PrimeFieldElement19 {
        PrimeFieldElement19::new(self.x.pow(pow).rem_euclid(MODULUS))
    }
}

impl Add<PrimeFieldElement19> for PrimeFieldElement19 {
    type Output = PrimeFieldElement19;

    fn add(self, other: PrimeFieldElement19) -> PrimeFieldElement19 {
        let x = (self.x + other.x).rem_euclid(MODULUS);
        PrimeFieldElement19::new(x)
    }
}

impl Sub<PrimeFieldElement19> for PrimeFieldElement19 {
    type Output = PrimeFieldElement19;

    fn sub(self, other: PrimeFieldElement19) -> PrimeFieldElement19 {
        let temp_x = self.x as i64;
        let temp_other = other.x as i64;
        let x: u64;
        if temp_other > temp_x {
            x = ((temp_x - temp_other + MODULUS as i64) as u64).rem_euclid(MODULUS);
        } else {
            x = ((temp_x - temp_other as i64) as u64).rem_euclid(MODULUS);
        }

        PrimeFieldElement19::new(x)
    }
}

impl Mul<PrimeFieldElement19> for PrimeFieldElement19 {
    type Output = PrimeFieldElement19;

    fn mul(self, other: PrimeFieldElement19) -> PrimeFieldElement19 {
        let x = (self.x * other.x).rem_euclid(MODULUS);
        PrimeFieldElement19::new(x)
    }
}

impl Mul<u64> for PrimeFieldElement19 {
    type Output = PrimeFieldElement19;

    fn mul(self, other: u64) -> PrimeFieldElement19 {
        let x = (self.x * other.rem_euclid(MODULUS)).rem_euclid(MODULUS);
        PrimeFieldElement19::new(x)
    }
}

impl Mul<PrimeFieldElement19> for u64 {
    type Output = PrimeFieldElement19;

    fn mul(self, other: PrimeFieldElement19) -> PrimeFieldElement19 {
        //self * other
        other * self
    }
}

// Inverses
impl Neg for PrimeFieldElement19 {
    type Output = PrimeFieldElement19;

    fn neg(self) -> PrimeFieldElement19 {
        PrimeFieldElement19::new(extended_gcd(self.x, MODULUS).rem_euclid(MODULUS))
    }
}

pub fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        return a;
    } else {
        return gcd(b, a.rem_euclid(b));
    }
}

/// Compute inverese of a modulo n using Extended Euclidean Algorithm.
pub fn extended_gcd(a: u64, n: u64) -> u64 {
    let mut t = 0i64;
    let mut r = n as i64;
    let mut newt = 1i64;
    let mut newr = a as i64;
    let mut quotient: i64;
    let mut temp: i64;

    while newr != 0 {
        quotient = r / newr;

        temp = newt;
        newt = t - quotient * newt;
        t = temp;

        temp = newr;
        newr = r - quotient * newr;
        r = temp;
    }

    if r > 1 {
        panic!("a is not invertible");
    } else if t < 0 {
        t = t + n as i64;
    }
    t as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd() {
        assert_eq!(3u64, gcd(9u64, 12u64));
        assert_eq!(1u64, gcd(5u64, 11u64));
    }

    #[test]
    fn test_extended_gcd() {
        assert_eq!(10u64, extended_gcd(2u64, 19u64));
    }

    #[test]
    fn test_inverse_of_field_element() {
        let p1 = PrimeFieldElement19::new(2);
        assert_eq!(PrimeFieldElement19::new(10), -p1);
    }

    #[test]
    fn add_field_elements() {
        let p1 = PrimeFieldElement19::new(5);
        let p2 = PrimeFieldElement19::new(5);
        let result = p1 + p2;
        assert_eq!(result.x, 10);
    }

    #[test]
    fn sub_field_elements() {
        let p1 = PrimeFieldElement19::new(10);
        let p2 = PrimeFieldElement19::new(20);
        let result = p1 - p2;
        assert_eq!(result.x, 9);

        let p1 = PrimeFieldElement19::new(16);
        let p2 = PrimeFieldElement19::new(30);
        let result = p1 - p2;
        assert_eq!(result.x, 5);

        let p1 = PrimeFieldElement19::new(5);
        let p2 = PrimeFieldElement19::new(15);
        let result = p1 - p2;
        assert_eq!(result.x, 9);
    }

    #[test]
    fn mul_field_elements() {
        let p1 = PrimeFieldElement19::new(1);
        let p2 = PrimeFieldElement19::new(20);
        let result = p1 * p2;
        assert_eq!(result.x, 1);
    }
}
