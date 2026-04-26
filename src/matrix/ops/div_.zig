const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const matops = @import("../ops.zig");

/// Performs in-place computation of the division of a matrix `x` and a numeric
/// `y` into a matrix `o`.
///
/// ## Signature
/// ```zig
/// matrix.div_(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the
///   same dimensions.
pub fn div_(o: anytype, x: anytype, y: anytype) !void {
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isMatrix(meta.Child(O)) or
        !meta.isMatrix(X) or !meta.isNumeric(Y))
        @compileError("zsl.matrix.div_: o must be a mutable one-itme pointer to a matrix, x must be a matrix, and y must be a numeric, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return matops.apply2_(o, x, y, numeric.div_);
}
