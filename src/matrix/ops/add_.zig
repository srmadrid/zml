const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const matops = @import("../ops.zig");

/// Performs in-place computation of the addition of two matrices `x` and `y`
/// into a matrix `o`.
///
/// ## Signature
/// ```zig
/// matrix.add_(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right matrix operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the three matrices do not have the
///   same dimensions.
pub fn add_(o: anytype, x: anytype, y: anytype) !void {
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isMatrix(meta.Child(O)) or
        !meta.isMatrix(X) or !meta.isMatrix(Y))
        @compileError("zsl.matrix.add_: o must be a mutable one-itme pointer to a matrix, and x and y must be matrices, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return matops.apply2_(o, x, y, numeric.add_);
}
