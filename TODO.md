# To Do

Maybe make it so shape and strides are stack-allocated. This way `allocator` can be optional and views (which do not need allocators as they access preexistent data) can be created without the need of an allocator and do not need to be deallocated.

Make two inits: `init` without flags (default values) and `initFlags` with flags.

Make `NDArray` agnostic to the type of the elements.

Make use of global context? Like, we have a type `ZMLContext` that has to be initialized before using the symbolic system using something like `zml.Context.init(allocator, other stuff)`. Then you have to set a global context, like `zml.SetGlobalContext(&context)`.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use setFromExpression(...). The expression must evaluate to the correct type.

Set.cartesian(null, &[_]*Set{&S, &S}), returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.
