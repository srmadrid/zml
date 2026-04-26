const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Min(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.min: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Min",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.min: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Min(type, type) type`");

            return Impl.Min(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Min", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.min: " ++ @typeName(X) ++ " must implement `fn Min(type, type) type`");

            return X.Min(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Min", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.min: " ++ @typeName(Y) ++ " must implement `fn Min(type, type) type`");

        return Y.Min(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return bool,
            .int => return int.Min(X, Y),
            .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .complex => @compileError("zsl.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Min(X, Y),
            .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .complex => @compileError("zsl.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .complex => @compileError("zsl.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Min(X, Y),
            .complex => @compileError("zsl.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zsl.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}

/// Returns the minimum between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.min(x: X, y: Y) numeric.Min(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Min(@TypeOf(x), @TypeOf(y))`: The minimum between `x` and `y`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Min` method. The expected
/// signature and behavior of `Min` are as follows:
/// * `fn Min(type, type) type`: Returns the type of `x + y`.
///
/// `numeric.Min(X, Y)`, `X` or `Y` must implement the required `min` method.
/// The expected signatures and behavior of `min` are as follows:
/// * `fn min(X, Y) numeric.Min(X, Y)`: Returns the minimum between `x` and
///   `y`.
pub fn min(x: anytype, y: anytype) numeric.Min(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Min(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "min",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.min: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.min(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "min",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.min(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "min",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.min(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return x and y,
            .int => return int.min(x, y),
            .float => return float.min(x, y),
            .dyadic => return dyadic.min(x, y),
            .complex => unreachable,
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.min(x, y),
            .float => return float.min(x, y),
            .dyadic => return dyadic.min(x, y),
            .complex => unreachable,
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.min(x, y),
            .dyadic => return dyadic.min(x, y),
            .complex => unreachable,
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.min(x, y),
            .complex => unreachable,
            .custom => unreachable,
        },
        .complex => unreachable,
        .custom => unreachable,
    }
}
