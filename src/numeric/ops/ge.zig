const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Compares any two numerics `x` and `y` for greater-than or equal ordering.
///
/// ## Signature
/// ```zig
/// numeric.ge(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than or equal to `y`, `false` otherwise.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `ge` method. The expected
/// signature and behavior of `ge` are as follows:
/// * `fn ge(X, Y) bool`: Compares `x` and `y` for greater-than or equal
///   ordering.
pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.ge: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "ge",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.ge: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn ge(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.ge(x, y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "ge", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zsl.numeric.ge: " ++ @typeName(X) ++ " must implement `fn ge(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.ge(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "ge", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zsl.numeric.ge: " ++ @typeName(Y) ++ " must implement `fn ge(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.ge(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return x or !y,
            .int => return int.ge(x, y),
            .float => return float.ge(x, y),
            .dyadic => return dyadic.ge(x, y),
            .complex => @compileError("zsl.numeric.ge: not defiged for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.ge(x, y),
            .float => return float.ge(x, y),
            .dyadic => return dyadic.ge(x, y),
            .complex => @compileError("zsl.numeric.ge: not defiged for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.ge(x, y),
            .dyadic => return dyadic.ge(x, y),
            .complex => @compileError("zsl.numeric.ge: not defiged for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.ge(x, y),
            .complex => @compileError("zsl.numeric.ge: not defiged for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zsl.numeric.ge: not defiged for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}
