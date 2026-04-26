const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Compares any two numerics `x` and `y` for less-than ordering.
///
/// ## Signature
/// ```zig
/// numeric.lt(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than `y`, `false` otherwise.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `lt` method. The expected
/// signature and behavior of `lt` are as follows:
/// * `fn lt(X, Y) bool`: Compares `x` and `y` for less-than ordering.
pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.lt: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "lt",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.lt: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.lt(x, y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "lt", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zsl.numeric.lt: " ++ @typeName(X) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.lt(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "lt", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zsl.numeric.lt: " ++ @typeName(Y) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.lt(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return !x and y,
            .int => return int.lt(x, y),
            .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .complex => @compileError("zsl.numeric.lt: not defiltd for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.lt(x, y),
            .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .complex => @compileError("zsl.numeric.lt: not defiltd for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .complex => @compileError("zsl.numeric.lt: not defiltd for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.lt(x, y),
            .complex => @compileError("zsl.numeric.lt: not defiltd for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zsl.numeric.lt: not defiltd for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}
