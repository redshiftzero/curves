//! Toy implementations of elliptic curves.
//!
//! ## References:
//!
//! * RFC 6090: Fundamental Elliptic Curve Cryptography Algorithms - `https://datatracker.ietf.org/doc/html/rfc6090`
//! * RFC 7748: Elliptic Curves for Security - `https://datatracker.ietf.org/doc/html/rfc7748`

pub mod edwards;
pub mod field;
pub mod scalar;
pub mod traits;
pub mod twisted_edwards;
pub mod weierstrass;
