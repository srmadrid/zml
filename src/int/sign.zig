const types = @import("../types.zig");

/// Returns the sign of an int `x`.
///
/// ## Signature
/// ```zig
/// int.sign(x: X) X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The int value to get the sign of.
///
/// ## Returns
/// `@TypeOf(x)`: The sign of `x`.
pub inline fn sign(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .int)
        @compileError("zml.int.sign: x must be an int, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (comptime @typeInfo(X).int.signedness == .unsigned) {
        return if (x == 0) 0 else 1;
    } else {
        return if (x > 0) 1 else if (x < 0) -1 else 0;
    }
}
