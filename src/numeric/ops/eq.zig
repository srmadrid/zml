const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Compares any two numerics `x` and `y` for equality.
///
/// ## Signature
/// ```zig
/// numeric.eq(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are equal, `false` otherwise.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `eq` method. The expected
/// signature and behavior of `eq` are as follows:
/// * `fn eq(X, Y) bool`: Compares `x` and `y` for equality.
pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.eq: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "eq",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.eq: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.eq(x, y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "eq", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zsl.numeric.eq: " ++ @typeName(X) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.eq(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "eq", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zsl.numeric.eq: " ++ @typeName(Y) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.eq(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return x == y,
            .int => return int.eq(x, y),
            .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .complex => switch (comptime meta.numericType(meta.Scalar(Y))) {
                .bool, .int => unreachable,
                .float => return float.eq(x, y.re) and float.eq(numeric.zero(X), y.im),
                .dyadic => return dyadic.eq(x, y.re) and dyadic.eq(numeric.zero(X), y.im),
                .complex, .custom => unreachable,
            },
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.eq(x, y),
            .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .complex => switch (comptime meta.numericType(meta.Scalar(Y))) {
                .bool, .int => unreachable,
                .float => return float.eq(x, y.re) and float.eq(numeric.zero(X), y.im),
                .dyadic => return dyadic.eq(x, y.re) and dyadic.eq(numeric.zero(X), y.im),
                .complex, .custom => unreachable,
            },
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .complex => switch (comptime meta.numericType(meta.Scalar(Y))) {
                .bool, .int => unreachable,
                .float => return float.eq(x, y.re) and float.eq(numeric.zero(X), y.im),
                .dyadic => return dyadic.eq(x, y.re) and dyadic.eq(numeric.zero(X), y.im),
                .complex, .custom => unreachable,
            },
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.eq(x, y),
            .complex => switch (comptime meta.numericType(meta.Scalar(Y))) {
                .bool, .int => unreachable,
                .float => return float.eq(x, y.re) and float.eq(numeric.zero(X), y.im),
                .dyadic => return dyadic.eq(x, y.re) and dyadic.eq(numeric.zero(X), y.im),
                .complex, .custom => unreachable,
            },
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(meta.Scalar(X))) {
            .bool, .int => unreachable,
            .float => switch (comptime meta.numericType(Y)) {
                .bool, .int, .float => return float.eq(x.re, y) and float.eq(x.im, numeric.zero(Y)),
                .dyadic => return dyadic.eq(x.re, y) and dyadic.eq(x.im, numeric.zero(Y)),
                .complex => switch (comptime meta.numericType(meta.Scalar(Y))) {
                    .bool, .int => unreachable,
                    .float => return float.eq(x.re, y.re) and float.eq(x.im, y.im),
                    .dyadic => return dyadic.eq(x.re, y.re) and dyadic.eq(x.im, y.im),
                    .complex, .custom => unreachable,
                },
                .custom => unreachable,
            },
            .dyadic => switch (comptime meta.numericType(Y)) {
                .bool, .int, .float, .dyadic => return dyadic.eq(x.re, y) and dyadic.eq(x.im, numeric.zero(Y)),
                .complex => return dyadic.eq(x.re, y.re) and dyadic.eq(x.im, y.im),
                .custom => unreachable,
            },
            .complex, .custom => unreachable,
        },
        .custom => unreachable,
    }
}
