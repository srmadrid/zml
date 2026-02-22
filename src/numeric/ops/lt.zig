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
/// ### `X` is custom
/// `X` must implement the required `lt` method. The expected signature and
/// behavior of `lt` are as follows:
/// * `fn lt(X, Y) bool`: Compares `x` and `y` for less-than ordering.
///
/// ### `Y` is a custom numeric type
/// `Y` must implement the required `lt` method, with the same specifications as
/// above.
///
/// ### Both `X` and `Y` are custom numeric types
/// At least one of `X` and `Y` must implement the required `lt` method, with
/// the same specifications as above. If both implement it, `X`'s implementation
/// will be used.
pub inline fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.haveMethod(
                &.{ X, Y },
                "lt",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.lt: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.lt(x, y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "lt", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zml.numeric.lt: " ++ @typeName(X) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.lt(x, y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "lt", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zml.numeric.lt: " ++ @typeName(Y) ++ " must implement `fn lt(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.lt(x, y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return !x and y,
            .int => return int.lt(x, y),
            .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.lt(x, y),
            .rational => return rational.lt(x, y),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.lt(x, y),
            .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.lt(x, y),
            .rational => return rational.lt(x, y),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.lt(x, y),
            .dyadic => return dyadic.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.lt(x, y),
            .rational => return rational.lt(x, y),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.lt(x, y),
            .rational => return rational.lt(x, y),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer => return integer.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .rational => return rational.lt(x, y),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer, .rational => return rational.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .real => return real.lt(x, y),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer, .rational, .real => return real.lt(x, y),
            .cfloat => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zml.numeric.lt: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}
