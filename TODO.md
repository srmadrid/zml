# To Do

## Priority

Complex `rotg` does not work. `rotmg` fails when compiling with optimization `ReleaseFast`.

Check all BLAS functions are correctly tested (especially the old ones: replace with the 2 and -2 incx tests). From `syr2` they follow the new (preferred) format and, therefore, have correct tests.

When all BLAS and LAPACK functions are implemented, clean `NDArray`: remove unnecessary functions, and correct any wrong ones (focus especially on the newly added `offset` parameter). Do the same for the iterators. Maybe remove `allocator` from `NDArray` and ask for it only when needed (currently only needed for `init` and `deinit`). For instance, for the arithmetic operations, use the template `op(allocator: ?Allocator, A: NDArray, B: NDArray) NDArray`: if `allocator` is null, then the operation is done in place on `A`, and a shallow copy is returned; if not, a new `NDArray` is created with the given allocator (if done in place, `A`'s dimensions must be appropriate). Maybe just have `op` and `opInPlace` functions. Maybe use pytorch's style: Methods which mutate a tensor are marked with an underscore suffix. For example, torch.FloatTensor.abs_() computes the absolute value in-place and returns the modified tensor, while torch.FloatTensor.abs() computes the result in a new tensor.

When changed, replace all `@setRuntimeSafety(false)` with `@optimizeFor(.ReleaseFast)`.

Eventually move the documentation website to a custom website (maybe using something like doxygen) instead of using zig's documentation system.

## General

Make a C interface (and, therefore, a C++ one) for the library.

Make a Python package (especially for the symbolic system).

Allow the user to build only BLAS and LAPACK to create `libblas.so`, `liblapack.so`, `libblas.a` and `liblapack.a` (also for Windows and MacOS) and export the necessary headers.

## `NDArray`

Add view functions: `reshape`, `slice`.

Offer storage options for `NDArray` (e.g., `Packed`, `Sparse` (`CSR`, `CSC`, `COO`), `Diagonal`, `Triangular`, `Symmetric`, `Hermitian`, etc.) instead of having different types for each one? If yes, a metadata field must be added to `NDArray` to store the extra information. Also, rename `RowMajor` and `ColumnMajor` to `RowMajorContiguous` and `ColumnMajorContiguous`, respectively.

Make two inits: `init` without flags (default values) and `initFlags` with flags?

Make `NDArray` agnostic to the type of the elements.

## Symbolic System

Make use of global context? Like, we have a type `zml.Context` that has to be initialized before using the symbolic system using something like `zml.Context.init(allocator, other stuff)`. Then you have to set a global context, like `zml.setGlobalContext(&context)`. The context will hold a `StringArrayHashMap` (hash map optimized for iteration) with symbols and their names. This way, the context stores all the symbols and the user can have several contexts at the same time and switch between them. The context will also store the allocator and other stuff.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use `setFromExpression(...)`. The expression must evaluate to the correct type.

`Set.cartesian(null, &[_]*Set{&S, &S})`, returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.
