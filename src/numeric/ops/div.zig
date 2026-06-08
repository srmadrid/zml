const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Div(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.Div: X and Y must be numeric types, got \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Div",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.Div: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Div(type, type) type`");

            return Impl.Div(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Div", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.Div: " ++ @typeName(X) ++ " must implement `fn Div(type, type) type`");

            return X.Div(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Div", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.Div: " ++ @typeName(Y) ++ " must implement `fn Div(type, type) type`");

        return Y.Div(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => @compileError("zsl.numeric.Div: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => return int.Div(X, Y),
            .float => return float.Div(X, Y),
            .dyadic => return dyadic.Div(X, Y),
            .complex => return complex.Div(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Div(X, Y),
            .float => return float.Div(X, Y),
            .dyadic => return dyadic.Div(X, Y),
            .complex => return complex.Div(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Div(X, Y),
            .dyadic => return dyadic.Div(X, Y),
            .complex => return complex.Div(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Div(X, Y),
            .complex => return complex.Div(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Div(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs division between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.div(x: X, y: Y) numeric.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Div` method. The expected
/// signature and behavior of `Div` are as follows:
/// * `fn Div(type, type) type`: Returns the type of `x/y`.
///
/// `numeric.Div(X, Y)`, `X` or `Y` must implement the required `div` method.
/// The expected signatures and behavior of `div` are as follows:
/// * `fn div(X, Y) numeric.Div(X, Y)`: Returns the division of `x` and `y`.
pub fn div(x: anytype, y: anytype) numeric.Div(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Div(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "div",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.div: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn div(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.div(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "div",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.div: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn div(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.div(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "div",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.div: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn div(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.div(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => unreachable,
            .int => return int.div(x, y),
            .float => return float.div(x, y),
            .dyadic => return dyadic.div(x, y),
            .complex => return complex.div(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.div(x, y),
            .float => return float.div(x, y),
            .dyadic => return dyadic.div(x, y),
            .complex => return complex.div(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.div(x, y),
            .dyadic => return dyadic.div(x, y),
            .complex => return complex.div(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.div(x, y),
            .complex => return complex.div(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.div(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs computation of the division of two numerics `x` and `y` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.divInto(o: *O, x: X, y: Y) void
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
/// `O`, `X` or `Y` should implement the required `divInto` method. The expected
/// signature and behavior of `divInto` are as follows:
/// * `fn divInto(*O, X, Y) void`: Computes the division of `x` and `y` and
///   stores it in `o`.
///
/// If none of `O`, `X` and `Y` implement the required `divInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.div`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X` and `Y` must adhere to the requirements of these
/// functions.
pub fn divInto(o: anytype, x: anytype, y: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y))
        @compileError("zsl.numeric.divInto: o must be a mutable one-item pointer to a numeric, and x and y must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // O, X and Y all custom
                if (comptime meta.anyHasMethod(&.{ O, X, Y }, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.divInto(o, x, y);
            } else { // only O and X custom
                if (comptime meta.anyHasMethod(&.{ O, X }, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.divInto(o, x, y);
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) { // only O and Y custom
                if (comptime meta.anyHasMethod(&.{ O, Y }, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.divInto(o, x, y);
            } else { // only O custom
                if (comptime meta.hasMethod(O, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return O.divInto(o, x, y);
            }
        }
    } else {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.divInto(o, x, y);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return X.divInto(o, x, y);
            }
        } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
            if (comptime meta.hasMethod(Y, "divInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                return Y.divInto(o, x, y);
        }
    }

    return numeric.set(o, numeric.div(x, y));
}
