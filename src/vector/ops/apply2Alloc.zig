const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

/// Applies a binary operation elementwise between two vectors, or between a
/// vector and a numeric, dynamically allocating memory for the result.
///
/// This function is intended for heap-allocated vector types (`vector.Dense`
/// and `vector.Sparse`), and performs dimension checks at runtime.
///
/// For two sparse vectors, or a sparse vector and a numeric, the operation is
/// only applied to the indices where at least one of the vectors has a non-zero
/// element.
///
/// ## Signature
/// ```zig
/// vector.apply2Allocated(allocator: std.mem.Allocator, x: X, y: Y, op: Op) !vector.Apply2(X, Y, op)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `op` (`comptime anytype`): A binary numeric function to apply elementwise
///   to `x` and `y`.
///
/// ## Returns
/// `vector.Apply2(@TypeOf(x), @TypeOf(y), op)`: The result of the operation.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `vector.Error.DimensionMismatch`: If the two vectors do not have the same
///   length. Can only happen if both operands are vectors.
pub fn apply2Alloc(allocator: std.mem.Allocator, x: anytype, y: anytype, comptime op: anytype) !vecops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Apply2(X, Y, op);

    const x_len_optional: ?usize = if (comptime meta.isVector(X)) (if (comptime meta.isStaticVector(X)) X.len else x.len) else null;
    const y_len = if (comptime meta.isVector(Y)) (if (comptime meta.isStaticVector(Y)) Y.len else y.len) else x_len_optional.?;
    const x_len = x_len_optional orelse y_len;

    if (comptime !((meta.isStaticVector(X) or meta.isNumeric(X)) and
        (meta.isStaticVector(Y) or meta.isNumeric(Y))))
    {
        if (x_len != y_len)
            return vector.Error.DimensionMismatch;
    }

    var result = switch (comptime meta.vectorType(R)) {
        .static => R.init,
        .dense => try R.init(allocator, x_len),
        .sparse => try R.init(allocator, x_len, (if (comptime meta.isSparseVector(X)) x.nnz else 0) + (if (comptime meta.isSparseVector(Y)) y.nnz else 0)),
        .numeric => unreachable,
    };

    vecops.apply2Into(
        &result,
        x,
        y,
        switch (comptime op) {
            numeric.add => numeric.addInto,
            numeric.sub => numeric.subInto,
            numeric.mul => numeric.mulInto,
            numeric.div => numeric.divInto,
            else => unreachable,
        },
    ) catch unreachable;

    return result;
}
