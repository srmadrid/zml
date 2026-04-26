const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Mul(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.mul: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Mul",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.mul: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Mul(type, type) type`");

            return Impl.Mul(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Mul", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.mul: " ++ @typeName(X) ++ " must implement `fn Mul(type, type) type`");

            return X.Mul(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Mul", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.mul: " ++ @typeName(Y) ++ " must implement `fn Mul(type, type) type`");

        return Y.Mul(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => @compileError("zsl.numeric.mul: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => return int.Mul(X, Y),
            .float => return float.Mul(X, Y),
            .dyadic => return dyadic.Mul(X, Y),
            .complex => return complex.Mul(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Mul(X, Y),
            .float => return float.Mul(X, Y),
            .dyadic => return dyadic.Mul(X, Y),
            .complex => return complex.Mul(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Mul(X, Y),
            .dyadic => return dyadic.Mul(X, Y),
            .complex => return complex.Mul(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Mul(X, Y),
            .complex => return complex.Mul(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Mul(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs multiplication between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.mul(x: X, y: Y) numeric.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Mul` method. The expected
/// signature and behavior of `Mul` are as follows:
/// * `fn Mul(type, type) type`: Returns the type of `x * y`.
///
/// `numeric.Mul(X, Y)`, `X` or `Y` must implement the required `mul` method.
/// The expected signatures and behavior of `mul` are as follows:
/// * `fn mul(X, Y) numeric.Mul(X, Y)`: Returns the multiplication of `x` and
///   `y`.
pub fn mul(x: anytype, y: anytype) numeric.Mul(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Mul(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "mul",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.mul: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn mul(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.mul(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "mul",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.mul: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn mul(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.mul(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "mul",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.mul: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn mul(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.mul(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => unreachable,
            .int => return int.mul(x, y),
            .float => return float.mul(x, y),
            .dyadic => return dyadic.mul(x, y),
            .complex => return complex.mul(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.mul(x, y),
            .float => return float.mul(x, y),
            .dyadic => return dyadic.mul(x, y),
            .complex => return complex.mul(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.mul(x, y),
            .dyadic => return dyadic.mul(x, y),
            .complex => return complex.mul(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.mul(x, y),
            .complex => return complex.mul(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.mul(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
