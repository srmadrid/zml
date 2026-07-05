const meta = @import("../meta.zig");

/// Returns the sign of an dyadic `x`.
///
/// ## Signature
/// ```zig
/// dyadic.sign(x: X) X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The dyadic value to get the sign of.
///
/// ## Returns
/// `@TypeOf(x)`: The sign of `x`.
pub fn sign(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or meta.numericType(X) != .dyadic)
        @compileError("zml.dyadic.sign: x must be an dyadic, got \n\tx: " ++ @typeName(X) ++ "\n");

    return if (x.isZero()) .zero else if (x.positive) .one else X.one.neg();
}
