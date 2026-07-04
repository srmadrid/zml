# zsl: zig scientific library

A generic numerical mathematics library for Zig.

> ⚠️ Zsl is in the early stages of development. APIs change frequently, results may be incorrect, and many features are missing or partial. Use with care, and please [open an issue](https://github.com/srmadrid/zsl/issues) if you hit a bug or have a suggestion.

## Why zsl?

Zsl lets you write numerical code once and run it over whatever numeric type your problem actually needs, e.g., `f32`, `f64`, dyadic rationals, complex numbers, or any custom numeric type, without incurring any unnecessary runtime cost.

## Current Features

- Numeric types:
  - Dyadic rationals
  - Complex numbers
- Vectors (`vector`):
  - Three storage formats:
    - `Static`
    - `Dense`
    - `Sparse`
  - Vector addition/subtraction and scalar multiplication/division
- Matrices (`matrix`):
  - Diverse storage formats:
    - General (`general`):
      - `Static`
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Symmetric (`symmetric`):
      - `Static`
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Hermitian (`hermitian`):
      - `Static`
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Triangular (`triangular`):
      - `Static`
      - `Dense`
      - `Sparse` (CSR, CSC)
    - Diagonal (`diagonal`):
      - `Static`
      - `Sparse`
    - Permutation (`permutation`):
      - `Static`
      - `Sparse`
  - Matrix addition/subtraction/multiplication and scalar multiplication/division
  - Views
- N-dimensional arrays (`array`):
  - Three storage formats:
    - `Static`
    - `Dense`
    - `Sparse` (CSF)
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
- Automatic Differentiation (`autodiff`):
  - Dual numbers (forward)
  - Tape and Var numbers (backward)

## Installation

To use this library in your project, run

```bash
zig fetch --save git+https://github.com/srmadrid/zsl
```

and add it to your build.zig file:

```zig
const zsl = b.dependency("zsl", .{});
exe.root_module.addImport("zsl", zsl.module("zsl"));
```

## Contributing
 
Contributions are very welcome. Feel free to open an issue or pr.
 
## License
 
zsl is released under the MIT License. See [LICENSE](LICENSE) for the full text.
