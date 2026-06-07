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
        @compileError("zsl.vector.Apply2: at least one of X or Y must be a vector type, the other must be a vector or a numeric type, and op must be a function of two arguments, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tOp = " ++ @typeName(Op) ++ "\n");

    comptime var R = meta.ReturnTypeFromInputs(op, &.{ meta.Numeric(X), meta.Numeric(Y) });
    const rinfo = @typeInfo(R);
    if (rinfo == .error_union)
        R = rinfo.error_union.payload;

    comptime if (!meta.isNumeric(R))
        @compileError("zsl.vector.Apply2: calling op with arguments of types " ++ @typeName(meta.Numeric(X)) ++ " and " ++ @typeName(meta.Numeric(Y)) ++ " must return a numeric, got\n\tR = " ++ @typeName(R) ++ "\n");

    switch (comptime meta.vectorType(X)) {
        .static => switch (comptime meta.vectorType(Y)) {
            .static => {
                if (comptime X.len != Y.len)
                    @compileError("zsl.vector.Apply2: static vector types X and Y must have equal lengths, got\n\tX = " ++
                        @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n");

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

/// Applies a binary operation elementwise between two static vectors, or
/// between a static vector and a numeric.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`), and performs dimension checks at compile time.
///
/// ## Signature
/// ```zig
/// vector.apply2(x: X, y: Y, op: Op) vector.Apply2(X, Y, op)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `op` (`comptime anytype`): A binary numeric function to apply elementwise
///   to `x` and `y`.
///
/// ## Returns
/// `vector.Apply2(@TypeOf(x), @TypeOf(y), op)`: The result of the operation.
pub fn apply2(x: anytype, y: anytype, comptime op: anytype) vecops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(op);
    const R: type = vecops.Apply2(X, Y, op);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X) or
        meta.isDenseVector(Y) or meta.isSparseVector(Y))
        @compileError("zsl.vector.apply2: the result cannot be a heap-allocated vector type, i.e., both inputs must be static vectors, or a static vector and a numeric, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n\tresult: " ++ @typeName(R) ++ "\n");

    var result = R.init;

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
