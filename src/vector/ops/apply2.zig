const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Apply2(comptime X: type, comptime Y: type, comptime op: anytype) type {
    const Op = @TypeOf(op);
    const opinfo = @typeInfo(Op);

    comptime if ((!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 2)
        @compileError("zsl.vector.apply2: at least one of x or y must be a vector, the other must be a vector or a numeric, and op must be a function of two arguments, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n");

    comptime var R = meta.ReturnTypeFromInputs(op, &.{ meta.Numeric(X), meta.Numeric(Y) });
    const rinfo = @typeInfo(R);
    if (rinfo == .error_union)
        R = rinfo.error_union.payload;

    comptime if (!meta.isNumeric(R))
        @compileError("zsl.vector.apply2: calling op with arguments of types X and Y must return a numeric, got\n\tR = " ++ @typeName(R) ++ "\n");

    switch (comptime meta.vectorType(X)) {
        .static => switch (comptime meta.vectorType(Y)) {
            .static => {
                if (comptime X.len != Y.len)
                    @compileError("zsl.vector.apply2: static vectors x and y must have equal compile time lengths, got\n\tx: " ++
                        @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n");

                return vector.Static(X.len, R);
            },
            .dense => return vector.Dense(R),
            .sparse => return vector.Dense(R),
            .numeric => return vector.Static(X.len, R),
        },
        .dense => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Dense(R),
            .dense => return vector.Dense(R),
            .sparse => return vector.Dense(R),
            .numeric => return vector.Dense(R),
        },
        .sparse => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Dense(R),
            .dense => return vector.Dense(R),
            .sparse => return vector.Sparse(R),
            .numeric => return vector.Sparse(R),
        },
        .numeric => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Static(Y.len, R),
            .dense => return vector.Dense(R),
            .sparse => return vector.Sparse(R),
            .numeric => unreachable,
        },
    }
}

/// Applies a binary operation elementwise between two vectors, or between a
/// vector and a numeric.
///
/// For two static vectors, the allocator is not used and can be set to
/// undefined, and the function cannot return an error unless op can.
///
/// For two sparse vectors, or a sparse vector and a numeric, the operation is
/// only applied to the indices where at least one of the vectors has a non-zero
/// element.
///
/// ## Signature
/// ```zig
/// vector.apply2(allocator: std.mem.Allocator, x: X, y: Y, op: Op) !vector.Apply2(X, Y, op)
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
pub fn apply2(allocator: std.mem.Allocator, x: anytype, y: anytype, comptime op: anytype) !vecops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Apply2(X, Y, op);

    const x_len = if (comptime meta.isVector(X)) x.len else y.len;
    const y_len = if (comptime meta.isVector(Y)) y.len else x.len;

    if (x_len != y_len)
        return vector.Error.DimensionMismatch;

    var result = switch (comptime meta.vectorType(R)) {
        .static => R.init,
        .dense => try R.init(allocator, x_len),
        .sparse => try R.init(allocator, x_len, (if (comptime meta.isSparseVector(X)) x.nnz else 0) + (if (comptime meta.isSparseVector(Y)) y.nnz else 0)),
        .numeric => unreachable,
    };

    vecops.apply2_(
        &result,
        x,
        y,
        switch (comptime op) {
            numeric.add => numeric.add_,
            numeric.sub => numeric.sub_,
            numeric.mul => numeric.mul_,
            numeric.div => numeric.div_,
            else => unreachable,
        },
    ) catch unreachable;

    return result;
}
