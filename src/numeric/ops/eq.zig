const types = @import("../../types.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const cfloat = @import("../../cfloat.zig");
const integer = @import("../../integer.zig");
const rational = @import("../../rational.zig");
const real = @import("../../real.zig");
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
/// ### `X` is custom
/// `X` must implement the required `eq` method. The expected signature and
/// behavior of `eq` are as follows:
/// * `fn eq(X, Y) bool`: Compares `x` and `y` for equality.
///
/// ### `Y` is a custom numeric type
/// `Y` must implement the required `eq` method, with the same specifications as
/// above.
///
/// ### Both `X` and `Y` are custom numeric types
/// At least one of `X` and `Y` must implement the required `eq` method, with
/// the same specifications as above. If both implement it, `X`'s implementation
/// will be used.
pub inline fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.haveMethod(
                &.{ X, Y },
                "eq",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.eq: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.eq(x, y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "eq", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zml.numeric.eq: " ++ @typeName(X) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.eq(x, y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "eq", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zml.numeric.eq: " ++ @typeName(Y) ++ " must implement `fn eq(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.eq(x, y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return x == y,
            .int => return int.eq(x, y),
            .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .cfloat => return cfloat.eq(x, y),
            .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.eq(x, y),
            .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .cfloat => return cfloat.eq(x, y),
            .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.eq(x, y),
            .dyadic => return dyadic.eq(x, y),
            .cfloat => return cfloat.eq(x, y),
            .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.eq(x, y),
            .cfloat => return cfloat.eq(x, y),
            .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat => return cfloat.eq(x, y),
            .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer => return integer.eq(x, y),
            .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational => return rational.eq(x, y),
            .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real => return real.eq(x, y),
            .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real, .complex => return complex.eq(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
