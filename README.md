# zml: zig mathematics library

A Zig numerical and symbolic mathematics library.

## Warning

This library is in the early stages of development and might return incorrect results. Breaking changes are to be expected every commit, and only the most basic functionality is currently implemented.

## Current Features

- Core:
  - Math:
    - Accurate math functions for all floating point types
  - Types:
    - Complex floats
    - Arbitrary precision integers, rationals, reals, and complex numbers: not implemented yet
- Vectors (`vector`):
  - Two storage formats:
    - `Dense`
    - `Sparse`
- Matrices (`matrix`):
  - Diverse storage formats:
    - General (`general`):
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Symmetric (`symmetric`):
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Hermitian (`hermitian`):
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Triangular (`triangular`):
      - `Dense`
      - `Sparse` (CSR, CSC)
    - `Diagonal`
    - `Permutation`
  - Matrix addition/subtraction and scalar multiplication/division
  - Views
- N-dimensional arrays (`array`):
  - Two storage formats:
    - `Dense` (plus `Strided` for views)
    - `Sparse` (CSF): not implemented yet
  - Broadcasting
  - Element-wise operations
  - Views
- Linear Algebra (`linalg`):
  - Matrix multiplication
  - Matrix decompositions:
    - LU (no pivoting (`lu`), partial pivoting (`plu`), full pivoting (`pluq`))
    - Cholesky (lower (`llt`), upper (`utu`), "smart" (`cholesky`))
    - Bunch-Kaufman (lower (`ldlt`), upper (`udut`), "smart" (`bunchkaufman`))
    - QR (no pivoting (`qr`), column pivoting (`qrp`))
  - BLAS routines (`blas`)
  - Select LAPACK routines (`lapack`)
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
