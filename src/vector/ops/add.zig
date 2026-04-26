const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.add: x and y must be vectors, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.add);
}

/// Performs addition between two vectors.
///
/// ## Signature
/// ```zig
/// vector.add(allocator: std.mem.Allocator, x: X, y: Y) !vector.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `vector.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `vector.Error.DimensionMismatch`: If the two vectors do not have the same
///   length.
pub fn add(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Add(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(allocator, x, y, numeric.add);
}
