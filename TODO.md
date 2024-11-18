# To Do

Make two inits: `init` without flags (default values) and `initFlags` with flags.

Make `NDArray` agnostic to the type of the elements.

Make a `<object>FromExpression` function. For instance, if I want to create the product of a set (`S\times S`), then I can use setFromExpression(...). The expression must evaluate to the correct type. 

Set.cartesian(null, &[_]*Set{&S, &S}), returns the set of cartesian product with "S \times S" as the default name (null was passed instead of a string). Same with other similar functions.
