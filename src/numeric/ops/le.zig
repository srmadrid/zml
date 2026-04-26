const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Compares any two numerics `x` and `y` for less-than or equal ordering.
///
/// ## Signature
/// ```zig
/// numeric.le(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than or equal to `y`, `false` otherwise.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `le` method. The expected
/// signature and behavior of `le` are as follows:
/// * `fn le(X, Y) bool`: Compares `x` and `y` for less-than or equal
///   ordering.
pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.le: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "le",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.le: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.le(x, y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "le", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zsl.numeric.le: " ++ @typeName(X) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.le(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "le", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zsl.numeric.le: " ++ @typeName(Y) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.le(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return !x or y,
            .int => return int.le(x, y),
            .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .complex => @compileError("zsl.numeric.le: not defiled for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.le(x, y),
            .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .complex => @compileError("zsl.numeric.le: not defiled for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .complex => @compileError("zsl.numeric.le: not defiled for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.le(x, y),
            .complex => @compileError("zsl.numeric.le: not defiled for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zsl.numeric.le: not defiled for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}
