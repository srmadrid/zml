# ZML: Zig Mathematics Library

A Zig numerical and symbolic mathematics library.

## Warning

This library is in the early stages of development and is not yet ready for use. Breaking changes are to be expected every commit, and only the most basic functionality is currently implemented.

## Current Features

- Core:
  - Math:
    - Accurate math functions for all floating point types
  - Types:
    - Complex floats
    - Arbitrary precision integers, rationals, reals, and complex numbers: not implemented yet
- Numerical System:
  - Matrices:
    - Diverse storage formats:
      - `General`
      - `Symmetric`
      - `Hermitian`
      - `Triangular`
      - `Diagonal`
      - `Banded`
      - `Tridiagonal`
      - `Sparse` (CSR, CSC): not implemented yet
    - Matrix addition/subtraction and scalar multiplication/division (in current development, only for select combinations)
    - Views
  - N-dimensional arrays:
    - Two storage formats:
      - `Dense` (plus `Strided` for views)
      - `Strided` (CSF)
    - Broadcasting
    - Element-wise operations
    - Views
  - Linear Algebra:
    - BLAS routines
    - Select LAPACK routines
- Symbolic System:
  - Nothing implemented yet

## Installation

To use this library in your project, run

```bash
zig fetch --save git+https://github.com/srmadrid/zml
```

and add it to your build.zig file:

```zig
const zml = b.dependency("zml", .{});
exe.root_module.addImport("zml", zml.module("zml"));
```

## Notes

This library is in the early stages of development and is not yet ready for use. Breaking changes are to be expected every commit, and only the most basic functionality is currently implemented (see [Current Features](#current-features) for more information). If you decide to use this library, please be aware that you may encounter bugs and incomplete features. If you find a bug, please report it via an issue on GitHub.

Any feature requests or suggestions are welcome. Please open an issue on GitHub to discuss any ideas you have.
