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
/// ### `X` is custom
/// `X` must implement the required `le` method. The expected signature and
/// behavior of `le` are as follows:
/// * `fn le(X, Y) bool`: Compares `x` and `y` for less-than or equal ordering.
///
/// ### `Y` is a custom numeric type
/// `Y` must implement the required `le` method, with the same specifications as
/// above.
///
/// ### Both `X` and `Y` are custom numeric types
/// At least one of `X` and `Y` must implement the required `le` method, with
/// the same specifications as above. If both implement it, `X`'s implementation
/// will be used.
pub inline fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.haveMethod(
                &.{ X, Y },
                "le",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.le: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.le(x, y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "le", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zml.numeric.le: " ++ @typeName(X) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.le(x, y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "le", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zml.numeric.le: " ++ @typeName(Y) ++ " must implement `fn le(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.le(x, y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return !x or y,
            .int => return int.le(x, y),
            .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.le(x, y),
            .rational => return rational.le(x, y),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.le(x, y),
            .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.le(x, y),
            .rational => return rational.le(x, y),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.le(x, y),
            .dyadic => return dyadic.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.le(x, y),
            .rational => return rational.le(x, y),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.le(x, y),
            .rational => return rational.le(x, y),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer => return integer.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .rational => return rational.le(x, y),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer, .rational => return rational.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .real => return real.le(x, y),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .integer, .rational, .real => return real.le(x, y),
            .cfloat => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zml.numeric.le: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}
