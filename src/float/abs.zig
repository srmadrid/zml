const types = @import("../types.zig");

/// Returns the absolute value of a float `x`.
///
/// ## Signature
/// ```zig
/// float.abs(x: X) X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The float value to get the absolute value of.
///
/// ## Returns
/// `@TypeOf(x)`: The absolute value of `x`.
pub inline fn abs(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.abs: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return @abs(x);
}
