# To Do

## Priority

Separate compile flags: `use_cblas` turn into `use_blas` (`bool`) to use BLAS functions in `NDArray` functions (inside use `options` to determine to use a BLAS function or a bespoke one), and, if true, have `link_cblas` (`[]const u8`) to link with a cblas library (if not provided, link with the internal zig implementation). Do the same with `use_clapack`.

Complex `rotg` does not work. `rotmg` fails when compiling with optimization `ReleaseFast`.

When changed, replace all `@setRuntimeSafety(false)` with `@optimizeFor(.ReleaseFast)`.

For `NDArray` add view functions: `reshape`, `flatten`, etc.

Edit `BLAS` functions parameters' in `NDArray` (pointers, const pointers, non pointers) to be consistent (keeping in mind that pinter or no pointer are just to show intent, as since we are editing the data buffer and not the struct itself, it does not matter if we pass a pointer or not). Using @constCast is a bad idea.

Maybe make strides `isize` instead of `usize` (probably not)?

## General

Make a C interface (and, therefore, a C++ one) for the library.

Make a Python package (especially for the symbolic system).

Allow the user to build only BLAS and LAPACK to create `libblas.so`, `liblapack.so`, `libblas.a` and `liblapack.a` (also for Windows and MacOS) and export the necessary headers.

Change the compilation option to use external BLAS and LAPACK libraries from bool to a string with the name of the library.

## `NDArray`

Offer two BLAS and LAPACK implementations: one for types that do not use allocators (bools, ints, floats and complex) and another for types that use allocators. The first one will be in the `BLAS` and `LAPACK` namespaces and the second one will be in the `BLASAlloc` and `LAPACKAlloc` namespaces (names not final). If the function signatures end up being the same, then use one namespace and use the `BLAS` and `LAPACK` names.

Maybe make it so shape and strides are stack-allocated. This way `allocator` can be optional and views (which do not need allocators as they access preexistent data) can be created without the need of an allocator and do not need to be deallocated.

Make two inits: `init` without flags (default values) and `initFlags` with flags?

Make `NDArray` agnostic to the type of the elements.

## Symbolic System

Make use of global context? Like, we have a type `zml.Context` that has to be initialized before using the symbolic system using something like `zml.Context.init(allocator, other stuff)`. Then you have to set a global context, like `zml.SetGlobalContext(&context)`. The context will hold a `StringArrayHashMap` (hash map optimized for iteration) with symbols and their names. This way, the context stores all the symbols and the user can have several contexts at the same time and switch between them. The context will also store the allocator and other stuff.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use setFromExpression(...). The expression must evaluate to the correct type.

Set.cartesian(null, &[_]*Set{&S, &S}), returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.
