const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.Add: X and Y must be vector types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.add);
}

/// Performs addition between two static vectors.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`), and performs dimension checks at compile time.
///
/// ## Signature
/// ```zig
/// vector.add(x: X, y: Y) vector.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `vector.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
pub fn add(x: anytype, y: anytype) vector.Add(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Add(X, Y);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X) or
        meta.isDenseVector(Y) or meta.isSparseVector(Y))
        @compileError("zsl.vector.add: the result cannot be a heap-allocated vector type, i.e., both inputs must be static vectors, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\n");

    return vecops.apply2(x, y, numeric.add);
}
