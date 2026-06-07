const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.Div: X must be a vector type and Y must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.div);
}

/// Performs division of a static vector by a numeric.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`).
///
/// ## Signature
/// ```zig
/// vector.div(x: X, y: Y) vector.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `vector.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) vector.Div(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Div(X, Y);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X))
        @compileError("zsl.vector.div: the result cannot be a heap-allocated vector type, i.e., x must be a static vector, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\n");

    return vecops.apply2(x, y, numeric.div);
}
