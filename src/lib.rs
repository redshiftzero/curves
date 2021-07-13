//! Toy implementations of elliptic curves.
//!
//! ## References:
//!
//! * RFC 6090: Fundamental Elliptic Curve Cryptography Algorithms - `https://datatracker.ietf.org/doc/html/rfc6090`
//! * RFC 7748: Elliptic Curves for Security - `https://datatracker.ietf.org/doc/html/rfc7748`

pub mod edwards;
pub mod field;
pub mod field_weierstrass;
pub mod real_weierstrass;
pub mod scalar;
pub mod traits;
