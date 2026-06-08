const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Hypot(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.Hypot: X and Y must be numeric types, got \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Hypot",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.Hypot: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Hypot(type, type) type`");

            return Impl.Hypot(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Hypot", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.Hypot: " ++ @typeName(X) ++ " must implement `fn Hypot(type, type) type`");

            return X.Hypot(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Hypot", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.Hypot: " ++ @typeName(Y) ++ " must implement `fn Hypot(type, type) type`");

        return Y.Hypot(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => @compileError("zsl.numeric.Hypot: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => @compileError("zsl.numeric.Hypot: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .float => return float.Hypot(X, Y),
            .dyadic => return dyadic.Hypot(X, Y),
            .complex => return complex.Hypot(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => @compileError("zsl.numeric.Hypot: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .float => return float.Hypot(X, Y),
            .dyadic => return dyadic.Hypot(X, Y),
            .complex => return complex.Hypot(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Hypot(X, Y),
            .dyadic => return dyadic.Hypot(X, Y),
            .complex => return complex.Hypot(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Hypot(X, Y),
            .complex => return complex.Hypot(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Hypot(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Computes the hypotenuse `√(x² + y²)` of any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.hypot(x: X, y: Y) numeric.Hypot(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Hypot(@TypeOf(x), @TypeOf(y))`: The hypotenuse of `x` and `y`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Hypot` method. The expected
/// signature and behavior of `Hypot` are as follows:
/// * `fn Hypot(type, type) type`: Returns the type of `√(x² + y²)`.
///
/// `numeric.Hypot(X, Y)`, `X` or `Y` must implement the required `hypot` method.
/// The expected signatures and behavior of `hypot` are as follows:
/// * `fn hypot(X, Y) numeric.Hypot(X, Y)`: Returns the hypotenuse of `x` and `y`.
pub fn hypot(x: anytype, y: anytype) numeric.Hypot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Hypot(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "hypot",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.hypot: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.hypot(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "hypot",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.hypot(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "hypot",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.hypot(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => unreachable,
            .int => unreachable,
            .float => return float.hypot(x, y),
            .dyadic => return dyadic.hypot(x, y),
            .complex => return complex.hypot(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => unreachable,
            .float => return float.hypot(x, y),
            .dyadic => return dyadic.hypot(x, y),
            .complex => return complex.hypot(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.hypot(x, y),
            .dyadic => return dyadic.hypot(x, y),
            .complex => return complex.hypot(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.hypot(x, y),
            .complex => return complex.hypot(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.hypot(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs computation of the hypotenuse `√(x² + y²)` of two numerics `x` and
/// `y` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.hypotInto(o: *O, x: X, y: Y) void
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
/// `O`, `X` or `Y` should implement the required `hypotInto` method. The
/// expected signature and behavior of `hypotInto` are as follows:
/// * `fn hypotInto(*O, X, Y) void`: Computes the hypotenuse of `x` and `y` and
///   stores it in `o`.
///
/// If none of `O`, `X` and `Y` implement the required `hypotInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.hypot`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X` and `Y` must adhere to the requirements of these
/// functions.
pub fn hypotInto(o: anytype, x: anytype, y: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y))
        @compileError("zsl.numeric.hypotInto: o must be a mutable one-item pointer to a numeric, and x and y must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // O, X and Y all custom
                if (comptime meta.anyHasMethod(&.{ O, X, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only O and X custom
                if (comptime meta.anyHasMethod(&.{ O, X }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) { // only O and Y custom
                if (comptime meta.anyHasMethod(&.{ O, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only O custom
                if (comptime meta.hasMethod(O, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return O.hypotInto(o, x, y);
            }
        }
    } else {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return X.hypotInto(o, x, y);
            }
        } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
            if (comptime meta.hasMethod(Y, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                return Y.hypotInto(o, x, y);
        }
    }

    return numeric.set(o, numeric.hypot(x, y));
}
