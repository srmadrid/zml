const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
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
/// `X` or `Y` must implement the required `ne` method. The expected
/// signature and behavior of `ne` are as follows:
/// * `fn ne(X, Y) bool`: Compares `x` and `y` for inequality.
pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.ne: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "ne",
                fn (X, Y) bool,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.ne: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return Impl.ne(x, y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "ne", fn (X, Y) bool, &.{ X, Y }))
                @compileError("zsl.numeric.ne: " ++ @typeName(X) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

            return X.ne(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "ne", fn (X, Y) bool, &.{ X, Y }))
            @compileError("zsl.numeric.ne: " ++ @typeName(Y) ++ " must implement `fn ne(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") bool`");

        return Y.ne(x, y);
    }

    return !numeric.eq(x, y);
}
