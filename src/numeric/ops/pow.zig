const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Pow(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.pow: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Pow",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.pow: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Pow(type, type) type`");

            return Impl.Pow(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Pow", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.pow: " ++ @typeName(X) ++ " must implement `fn Pow(type, type) type`");

            return X.Pow(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Pow", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.pow: " ++ @typeName(Y) ++ " must implement `fn Pow(type, type) type`");

        return Y.Pow(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => @compileError("zsl.numeric.pow: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => return int.Pow(X, Y),
            .float => return float.Pow(X, Y),
            .dyadic => return dyadic.Pow(X, Y),
            .complex => return complex.Pow(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Pow(X, Y),
            .float => return float.Pow(X, Y),
            .dyadic => return dyadic.Pow(X, Y),
            .complex => return complex.Pow(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Pow(X, Y),
            .dyadic => return dyadic.Pow(X, Y),
            .complex => return complex.Pow(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Pow(X, Y),
            .complex => return complex.Pow(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Pow(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs exponentiation `xʸ` between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.pow(x: X, y: Y) numeric.Pow(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Pow(@TypeOf(x), @TypeOf(y))`: The result of the exponentiation.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Pow` method. The expected
/// signature and behavior of `Pow` are as follows:
/// * `fn Pow(type, type) type`: Returns the type of `x + y`.
///
/// `numeric.Pow(X, Y)`, `X` or `Y` must implement the required `pow` method.
/// The expected signatures and behavior of `pow` are as follows:
/// * `fn pow(X, Y) numeric.Pow(X, Y)`: Returns the exponentiation of `x` to
///   the power `y`.
pub fn pow(x: anytype, y: anytype) numeric.Pow(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Pow(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "pow",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.pow: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn pow(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.pow(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "pow",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.pow: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn pow(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.pow(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "pow",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.pow: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn pow(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.pow(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => unreachable,
            .int => return int.pow(x, y),
            .float => return float.pow(x, y),
            .dyadic => return dyadic.pow(x, y),
            .complex => return complex.pow(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.pow(x, y),
            .float => return float.pow(x, y),
            .dyadic => return dyadic.pow(x, y),
            .complex => return complex.pow(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.pow(x, y),
            .dyadic => return dyadic.pow(x, y),
            .complex => return complex.pow(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.pow(x, y),
            .complex => return complex.pow(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.pow(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs computation of the exponentiation `xʸ` of two numerics `x` and `y`
/// into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.powInto(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O`, `X` or `Y` should implement the required `powInto` method. The expected
/// signature and behavior of `powInto` are as follows:
/// * `fn powInto(*O, X, Y) void`: Computes the exponentiation of `x` to the
///   power `y` and stores it in `o`.
///
/// If none of `O`, `X` and `Y` implement the required `powInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.pow`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X` and `Y` must adhere to the requirements of these
/// functions.
pub fn powInto(o: anytype, x: anytype, y: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y))
        @compileError("zsl.numeric.powInto: o must be a mutable one-item pointer to a numeric, and x and y must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // O, X and Y all custom
                if (comptime meta.anyHasMethod(&.{ O, X, Y }, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.powInto(o, x, y);
            } else { // only O and X custom
                if (comptime meta.anyHasMethod(&.{ O, X }, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.powInto(o, x, y);
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) { // only O and Y custom
                if (comptime meta.anyHasMethod(&.{ O, Y }, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.powInto(o, x, y);
            } else { // only O custom
                if (comptime meta.hasMethod(O, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return O.powInto(o, x, y);
            }
        }
    } else {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.powInto(o, x, y);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return X.powInto(o, x, y);
            }
        } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
            if (comptime meta.hasMethod(Y, "powInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                return Y.powInto(o, x, y);
        }
    }

    return numeric.set(o, numeric.pow(x, y));
}
