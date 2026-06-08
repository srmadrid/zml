const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Fma(X: type, Y: type, Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z))
        @compileError("zsl.numeric.Fma: X, Y and Z must be numeric types, got \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(Z)) { // X, Y and Z all custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ X, Y, Z },
                    "Fma",
                    fn (type, type, type) type,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.Fma: " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ " or " ++ @typeName(Z) ++ " must implement `fn Fma(type, type, type) type`");

                return Impl.Fma(X, Y, Z);
            } else { // only X and Y custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ X, Y },
                    "Fma",
                    fn (type, type, type) type,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.Fma: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Fma(type, type, type) type`");

                return Impl.Fma(X, Y, Z);
            }
        } else {
            if (comptime meta.isCustomNumeric(Z)) { // only X and Z custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ X, Z },
                    "Fma",
                    fn (type, type, type) type,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.Fma: " ++ @typeName(X) ++ " or " ++ @typeName(Z) ++ " must implement `fn Fma(type, type, type) type`");

                return Impl.Fma(X, Y, Z);
            } else { // only X custom
                comptime if (!meta.hasMethod(X, "Fma", fn (type, type, type) type, &.{ X, Y, Z }))
                    @compileError("zsl.numeric.Fma: " ++ @typeName(X) ++ " must implement `fn Fma(type, type, type) type`");

                return X.Fma(X, Y, Z);
            }
        }
    } else if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(Z)) { // only Y and Z custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ Y, Z },
                "Fma",
                fn (type, type, type) type,
                &.{ X, Y, Z },
            ) orelse
                @compileError("zsl.numeric.Fma: " ++ @typeName(Y) ++ " or " ++ @typeName(Z) ++ " must implement `fn Fma(type, type, type) type`");

            return Impl.Fma(X, Y, Z);
        } else { // only Y custom
            comptime if (!meta.hasMethod(Y, "Fma", fn (type, type, type) type, &.{ X, Y, Z }))
                @compileError("zsl.numeric.Fma: " ++ @typeName(Y) ++ " must implement `fn Fma(type, type, type) type`");

            return Y.Fma(X, Y, Z);
        }
    } else if (comptime meta.isCustomNumeric(Z)) { // only Z custom
        comptime if (!meta.hasMethod(Z, "Fma", fn (type, type, type) type, &.{ X, Y, Z }))
            @compileError("zsl.numeric.Fma: " ++ @typeName(Z) ++ " must implement `fn Fma(type, type, type) type`");

        return Z.Fma(X, Y, Z);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => switch (comptime meta.numericType(Z)) {
                .bool => @compileError("zsl.numeric.Fma: not defined for " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ " and " ++ @typeName(Z) ++ "."),
                .int => return int.Fma(X, Y, Z),
                .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .int => switch (comptime meta.numericType(Z)) {
                .bool, .int => return int.Fma(X, Y, Z),
                .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => switch (comptime meta.numericType(Z)) {
                .bool, .int => return int.Fma(X, Y, Z),
                .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.Fma(X, Y, Z),
                .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.Fma(X, Y, Z),
                .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.Fma(X, Y, Z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs fused multiplication and addition (x * y + z) between any three
/// numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.fma(x: X, y: Y, z: Z) numeric.Fma(X, Y, Z)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `numeric.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z))`: The result of the fused
/// multiplication and addition.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X`, `Y` or `Z` must implement the required `Fma` method. The expected
/// signature and behavior of `Fma` are as follows:
/// * `fn Fma(type, type, type) type`: Returns the type of `x * y + z`.
///
/// `numeric.Fma(X, Y, Z)`, `X`, `Y` or `Z` must implement the required `fma`
/// method. The expected signatures and behavior of `fma` are as follows:
/// * `fn fma(X, Y, Z) numeric.Fma(X, Y, Z)`: Returns the fused multiplication
/// and addition of `x`, `y` and `z`.
pub fn fma(x: anytype, y: anytype, z: anytype) numeric.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);
    const R: type = numeric.Fma(X, Y, Z);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(Z)) { // X, Y and Z all custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ R, X, Y, Z },
                    "fma",
                    fn (X, Y, Z) R,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ " or " ++ @typeName(Z) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

                return Impl.fma(x, y, z);
            } else { // only X and Y custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ R, X, Y },
                    "fma",
                    fn (X, Y, Z) R,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

                return Impl.fma(x, y, z);
            }
        } else {
            if (comptime meta.isCustomNumeric(Z)) { // only X and Z custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ R, X, Z },
                    "fma",
                    fn (X, Y, Z) R,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Z) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

                return Impl.fma(x, y, z);
            } else { // only X custom
                const Impl: type = comptime meta.anyHasMethod(
                    &.{ R, X },
                    "fma",
                    fn (X, Y, Z) R,
                    &.{ X, Y, Z },
                ) orelse
                    @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

                return Impl.fma(x, y, z);
            }
        }
    } else if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(Z)) { // only Y and Z custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, Y, Z },
                "fma",
                fn (X, Y, Z) R,
                &.{ X, Y, Z },
            ) orelse
                @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ ", " ++ @typeName(Y) ++ " or " ++ @typeName(Z) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.fma(x, y, z);
        } else { // only Y custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, Y },
                "fma",
                fn (X, Y, Z) R,
                &.{ X, Y, Z },
            ) orelse
                @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.fma(x, y, z);
        }
    } else if (comptime meta.isCustomNumeric(Z)) { // only Z custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Z },
            "fma",
            fn (X, Y, Z) R,
            &.{ X, Y, Z },
        ) orelse
            @compileError("zsl.numeric.fma: " ++ @typeName(R) ++ " or " ++ @typeName(Z) ++ " must implement `fn fma(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ", " ++ @typeName(Z) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.fma(x, y, z);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => switch (comptime meta.numericType(Z)) {
                .bool => unreachable,
                .int => return int.fma(x, y, z),
                .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .int => switch (comptime meta.numericType(Z)) {
                .bool, .int => return int.fma(x, y, z),
                .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => switch (comptime meta.numericType(Z)) {
                .bool, .int => return int.fma(x, y, z),
                .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float => return float.fma(x, y, z),
                .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic => return dyadic.fma(x, y, z),
                .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => switch (comptime meta.numericType(Z)) {
                .bool, .int, .float, .dyadic, .complex => return complex.fma(x, y, z),
                .custom => unreachable,
            },
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs computation of the fused multiplication and addition of three
/// numerics `x`, `y` and `z` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.fmaInto(o: *O, x: X, y: Y, z: Z) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O`, `X`, `Y` or `Z` should implement the required `fmaInto` method. The
/// expected signature and behavior of `fmaInto` are as follows:
/// * `fn fmaInto(*O, X, Y, Z) void`: Computes the fused multiplication and
///   addition of `x`, `y` and `z` and stores it in `o`.
///
/// If none of `O`, `X`, `Y` and `Z` implement the required `fmaInto` method,
/// the function will fall back to using `numeric.set` with the result of
/// `numeric.fma`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X`, `Y` and `Z` must adhere to the requirements of these
/// functions.
pub fn fmaInto(o: anytype, x: anytype, y: anytype, z: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y) or
        !meta.isNumeric(Z))
        @compileError("zsl.numeric.fmaInto: o must be a mutable one-item pointer to a numeric, and x, y and < must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tz: " ++ @typeName(Z) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) {
                if (comptime meta.isCustomNumeric(Z)) { // O, X, Y and Z all custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Y, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                } else { // only O, X and Y custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Y }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                }
            } else {
                if (comptime meta.isCustomNumeric(Z)) { // only O, X and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                } else { // only O and X custom
                    if (comptime meta.anyHasMethod(&.{ O, X }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                }
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) {
                if (comptime meta.isCustomNumeric(Z)) { // only O, Y and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, Y, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                } else { // only O and Y custom
                    if (comptime meta.anyHasMethod(&.{ O, Y }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                }
            } else {
                if (comptime meta.isCustomNumeric(Z)) { // only O and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fmaInto(o, x, y, z);
                } else { // only O custom
                    if (comptime meta.hasMethod(O, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                        return O.fmaInto(o, x, y, z);
                }
            }
        }
    } else if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(Z)) { // only X, Y and Z custom
                if (comptime meta.anyHasMethod(&.{ X, Y, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fmaInto(o, x, y, z);
            } else { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fmaInto(o, x, y, z);
            }
        } else {
            if (comptime meta.isCustomNumeric(Z)) { // only X and Z custom
                if (comptime meta.anyHasMethod(&.{ X, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fmaInto(o, x, y, z);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                    return X.fmaInto(o, x, y, z);
            }
        }
    } else if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(Z)) { // only Y and Z custom
            if (comptime meta.anyHasMethod(&.{ Y, Z }, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                return Impl.fmaInto(o, x, y, z);
        } else { // only Y custom
            if (comptime meta.hasMethod(Y, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                return Y.fmaInto(o, x, y, z);
        }
    } else if (comptime meta.isCustomNumeric(Z)) { // only Z custom
        if (comptime meta.hasMethod(Z, "fmaInto", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
            return Z.fmaInto(o, x, y, z);
    }

    return numeric.set(o, numeric.fma(x, y, z));
}
