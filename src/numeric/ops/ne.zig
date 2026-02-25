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

/// Compares any two numerics `x` and `y` for inequality.
///
/// ## Signature
/// ```zig
/// numeric.ne(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are not equal, `false` otherwise.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// ### `X` is custom
/// `X` must implement the required `ne` method. The expected signature and
/// behavior of `ne` are as follows:
/// * `fn ne(X, Y) bool`: Compares `x` and `y` for inequality.
///
/// ### `Y` is a custom numeric type
/// `Y` must implement the required `ne` method, with the same specifications as
/// above.
///
/// ### Both `X` and `Y` are custom numeric types
/// At least one of `X` and `Y` must implement the required `ne` method, with
/// the same specifications as above. If both implement it, `X`'s implementation
/// will be used.
pub inline fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.anyHasMethod(
                &.{ X, Y },
                "ne",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.ne: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.ne(x, y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "ne", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zml.numeric.ne: " ++ @typeName(X) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.ne(x, y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "ne", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zml.numeric.ne: " ++ @typeName(Y) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.ne(x, y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return x != y,
            .int => return int.ne(x, y),
            .float => return float.ne(x, y),
            .dyadic => return dyadic.ne(x, y),
            .cfloat => return cfloat.ne(x, y),
            .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.ne(x, y),
            .float => return float.ne(x, y),
            .dyadic => return dyadic.ne(x, y),
            .cfloat => return cfloat.ne(x, y),
            .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.ne(x, y),
            .dyadic => return dyadic.ne(x, y),
            .cfloat => return cfloat.ne(x, y),
            .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.ne(x, y),
            .cfloat => return cfloat.ne(x, y),
            .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat => return cfloat.ne(x, y),
            .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer => return integer.ne(x, y),
            .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational => return rational.ne(x, y),
            .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real => return real.ne(x, y),
            .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real, .complex => return complex.ne(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
