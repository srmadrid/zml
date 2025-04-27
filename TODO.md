# To Do

## Priority

Add inplace variants for core.math functions (e.g., `add_`, `mul_`, etc.) that take a pointer to the result, and check that the input types can be coerced to the output type.

Make exhaustive (test all input combinations) tests for type functions (core/types)

Implement basic math functions based on `ref/glibc/sysdeps/ieee754/ldbl-128ibm` IBM implementations.

Make coercion functions coerce ints into suitably big floats, and same for floats to cfloats.

Namespacing:

- `zml.add`: given any two objects of any type supported in the library (can be `NDArray` (or slice), `Set`, `Function`, etc.), executes it promoting the types to the most general type
- `zml.core.types.add`: given two numeric types, executes it promoting the types to the most general type
- `zml.core.types.cfloat`: given two fixed precision type variables (at least one must be cfloat), returns the most general cfloat.

This is the general idea

Create arbitrary precision types: Integer, Rational, Real and Complex.

Change `get` to return a pointer to the element instead of the element itself.

Since we only accept certain types (for both the symbolic and numerical systems), in core have the `add`, `sub`, `mul`, `div`, `neg`, `abs`, `pow`, `sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, etc. functions that add any two types and return the correct type, when possible. For instance, `add` can be used to add two `int`s, two `float`s, an `int` and a `float`, etc. The same for the other functions. This way, we can use them in the symbolic system and in the numerical system. In core.ops.

Create in ndarray/ops.zig simple operations. This will be used for functions like add(a: anytype, b: anytype) NDArray(Coerce(@TypeOf(a), @TypeOf(b))) {}, where Coerce(comptime T: anytype, comptime U: anytype) type returns the scalar type which can coerce from addition (and T and U can be NDArrays or scañars themselves).

Complex `rotg` does not work. `rotmg` fails when compiling with optimization `ReleaseFast`.

When all BLAS and LAPACK functions are implemented, clean `NDArray`: remove unnecessary functions, and correct any wrong ones (focus especially on the newly added `offset` parameter). Do the same for the iterators. Maybe remove `allocator` from `NDArray` and ask for it only when needed (currently only needed for `init` and `deinit` as there is never a need to realloc the array). For instance, for the arithmetic operations, use the template `op(allocator: ?Allocator, A: NDArray, B: NDArray) NDArray`: if `allocator` is null, then the operation is done in place on `A`, and a shallow copy is returned; if not, a new `NDArray` is created with the given allocator (if done in place, `A`'s dimensions must be appropriate). Maybe just have `op` and `opInPlace` functions. Maybe use pytorch's style: Methods which mutate a tensor are marked with an underscore suffix. For example, `torch.FloatTensor.abs_()` computes the absolute value in-place and returns the modified tensor, while `torch.FloatTensor.abs()` computes the result in a new tensor.

When changed, replace all `@setRuntimeSafety(false)` with `@optimizeFor(.ReleaseFast)`.

Eventually move the documentation website to a custom website (maybe using something like doxygen or docusaurus) instead of using zig's documentation system.

## General

Make a C interface (and, therefore, a C++ one) for the library.

Make a Python package (especially for the symbolic system).

Allow the user to build only BLAS and LAPACK to create `libblas.so`, `liblapack.so`, `libblas.a` and `liblapack.a` (also for Windows and MacOS) and export the necessary headers.

## `NDArray`

For functions which require optional parameters, use a `options` argument with a custom struct, like:

```rs
// Add won't actually take options, but this is just an example
pub fn add(
    allocator: Allocator,
    left: anytype,
    right: anytype,
    options: struct {
        alpha: f64 = 1.0,
        beta: f64 = 0.0,
        offsetA: usize = 0,
        offsetB: usize = 0,
        offsetC: usize = 0,
    },
) !NDArray(Coerce(Numeric(@TypeOf(left)), Numeric(@TypeOf(right)))) {
    // ...
}

// In place
pub fn add_(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        alpha: f64 = 1.0,
        beta: f64 = 0.0,
        offsetA: usize = 0,
        offsetB: usize = 0,
        offsetC: usize = 0,
    },
) !void {
    // ...
}
```

For the final all types included general function `zml.add`:

```rs
pub fn add(
    allocator: Allocator,
    left: anytype,
    right: anytype,
    options: struct {
        alpha: f64 = 1.0,
        beta: f64 = 0.0,
        offsetA: usize = 0,
        offsetB: usize = 0,
        offsetC: usize = 0,
    },
) !Coerce(@TypeOf(left), @TypeOf(right)) {
    // ...
}
```

In inplace functions, `out` should hold a type such that `canCoerce(@TypeOf(left), @TypeOf(right), @TypeOf(out))` is true.

All functions should take scalars, slices or NDArrays. If a slice is passed, call fromSlice and use that NDArray (assume RowMajor Packed)

`get` should return a pointer.

Add view functions: `reshape`, `slice`.

Since more storage options are being added, the `next()` function in iterators should take `axis` instead of `order`?

Offer storage options for `NDArray` (e.g., `Packed`, `Sparse` (`CSR`, `CSC`, `COO`), `Diagonal`, `Triangular`, `Symmetric`, `Hermitian`, etc.) instead of having different types for each one? If yes, a metadata field must be added to `NDArray` to store the extra information. Also, rename `RowMajor` and `ColumnMajor` to `RowMajorContiguous` and `ColumnMajorContiguous`, respectively. Maybe have `Order` with `RowMajor`, `ColumnMajor` and `Other` (not final name), and then `Storage` with `Packed` (`Packed` with row major or column major is like contiguous), `Sparse`, etc.? Also, Sparse Strided, when the data is not contiguous, but the indices are, so an efficient loop can still be made (like when using subarrays).

Make two inits: `init` without flags (default values) and `initFlags` with flags?

Test SPFM Lab7, Q1. Put as an example. Also maybe look at exercises from Statistical Signal Processing.

`initFn` will initialize all values with a function, like `initFn(...,normal.pdf(0,1),...)` creates an array with values extracted from a normal distribution.

## Symbolic System

Make use of global context? Like, we have a type `zml.Context` that has to be initialized before using the symbolic system using something like `zml.Context.init(allocator, other stuff)`. Then you have to set a global context, like `zml.setGlobalContext(&context)`. The context will hold a `StringArrayHashMap` (hash map optimized for iteration) with symbols and their names. This way, the context stores all the symbols and the user can have several contexts at the same time and switch between them. The context will also store the allocator and other stuff.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use `setFromExpression(...)`. The expression must evaluate to the correct type.

`Set.cartesian(null, &[_]*Set{&S, &S})`, returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.

For integrals, derivatives, etc. have anytype parameter inputs for the generic functions, and if the input is an expression of symbolic `Function` object, then calculate it symbolically, if it is a typical zig function, then calculate it numerically.

Have custom symbols, such as `\variance`, `\expectation`, etc., with some default values. For instance, `\variance` can be defined as `\newcommand{\variance}{\mathbb{V}\text{ar}}`, but that can be changed by the user to their liking.
