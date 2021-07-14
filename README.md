# curves

This repository contains various elliptic curve related examples. The implementations are not side-channel resistant, are not performant, and should not be used for anything other than educational purposes.

## Contents

### Affine space computations

* [Edwards curve over a small prime field (Rust)](src/edwards.rs)
* [Weierstrass curve over a small prime field (Rust)](src/weierstrass.rs)
* [Weierstrass curve over the reals (Python)](curves.py)

## Setup

Rust:

```
cargo test
cargo doc --open
```

Python:

```
python3 -m venv env
source env/bin/activate
python3 -m pip install -r requirements.txt
python3 curves.py
```
