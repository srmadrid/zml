# To Do

## Priority

For `NDArray` add view functions: `reshape`, `flatten`, `ravel`, `squeeze`, etc.

Edit `BLAS` functions parameters' (pointers, const pointers, non pointers) to be consistent (keeping in mind that pinter or no pointer are just to show intent, as since we are editing the data buffer and not the struct itself, it does not matter if we pass a pointer or not). Using @constCast is a bad idea.

Make BLAS functions take `T` and `T`, and inside `NDArray(T)` just make calls to them. This way all BLAS and LAPACK functions are according to the specification on netlib.

Maybe make strides `isize` instead of `usize`?

## General

Make a C interface.

Somehow be able to apply some optimization to some functions, and a different optimization to the rest (basically, force `BLAS` and `LAPACK` to be ReleaseFast, and the rest to be whatever the user wants).

## `NDArray`

Maybe, for vector and matrix BLAS functions, expect correct shapes, and the user can just do a `reshape` if needed (`flatten` is a view, so it does not need to be deallocated, but needs to be able to "shrink" the shape).

In the build system, give the option for the user to choose the BLAS and LAPACK implementation. If the user does not choose, then use the default one implemented in the library using zig.

Offer two BLAS and LAPACK implementations: one for types that do not use allocators (bools, ints, floats and complex) and another for types that use allocators. The first one will be in the `BLAS` and `LAPACK` namespaces and the second one will be in the `BLASAlloc` and `LAPACKAlloc` namespaces (names not final). If the function signatures end up being the same, then use one namespace and use the `BLAS` and `LAPACK` names.

Maybe make it so shape and strides are stack-allocated. This way `allocator` can be optional and views (which do not need allocators as they access preexistent data) can be created without the need of an allocator and do not need to be deallocated.

Make two inits: `init` without flags (default values) and `initFlags` with flags?

Make `NDArray` agnostic to the type of the elements.

## Symbolic System

Make use of global context? Like, we have a type `zml.Context` that has to be initialized before using the symbolic system using something like `zml.Context.init(allocator, other stuff)`. Then you have to set a global context, like `zml.SetGlobalContext(&context)`. The context will hold a `StringArrayHashMap` with symbols and their names. This way, the context stores all the symbols and the user can have several contexts at the same time and switch between them. The context will also store the allocator and other stuff.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use setFromExpression(...). The expression must evaluate to the correct type.

Set.cartesian(null, &[_]*Set{&S, &S}), returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.
