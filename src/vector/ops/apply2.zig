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

    comptime if (op != numeric.add and op != numeric.sub and op != numeric.mul and op != numeric.div)
        @compileError("zsl.vector.Apply2: op must be zsl.numeric.add, zsl.numeric.sub, zsl.numeric.mul or zsl.numeric.div, got\n\top: " ++ @typeName(Op) ++ "\n");

    comptime if (meta.isStaticVector(X) and meta.isStaticVector(Y) and X.len != Y.len)
        @compileError("zsl.vector.Apply2: static vector types X and Y must have equal lengths, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n");

    switch (comptime meta.vectorType(X)) {
        .static => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Static(X.len, R),
            .dense => return vector.Static(X.len, R),
            .sparse => return vector.Static(X.len, R),
            .numeric => return vector.Static(X.len, R),
        },
        .dense => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Static(Y.len, R),
            .dense => return vector.Dense(R),
            .sparse => return vector.Dense(R),
            .numeric => return vector.Dense(R),
        },
        .sparse => switch (comptime meta.vectorType(Y)) {
            .static => return vector.Static(Y.len, R),
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
/// This function is intended for when the result vector's length is known at
/// compile time, i.e., at least one of the inputs is a static vector. For two
/// static vectors, or a static vector and a numeric dimension checks are
/// performed at compile time, for any other combination, dimension checks are
/// performed at runtime throught `std.debug.assert`.
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

    if (comptime meta.isDenseVector(R) or meta.isSparseVector(R))
        @compileError("zsl.vector.apply2: the result cannot be a heap-allocated vector type, i.e., at least one input must be a static vector, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.vector.apply2Alloc instead.");

    const x_len_optional: ?usize = if (comptime meta.isVector(X)) (if (comptime meta.isStaticVector(X)) X.len else x.len) else null;
    const y_len = if (comptime meta.isVector(Y)) (if (comptime meta.isStaticVector(Y)) Y.len else y.len) else x_len_optional.?;
    const x_len = x_len_optional orelse y_len;

    if (comptime !(meta.isStaticVector(X) or meta.isNumeric(X)) or !(meta.isStaticVector(Y) or meta.isNumeric(Y)))
        std.debug.assert(x_len == y_len);

    var result = R.init;

    vecops.apply2IntoUnchecked(
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
    );

    return result;
}

/// Applies a binary operation elementwise between two vectors, or between a
/// vector and a numeric, dynamically allocating memory for the result.
///
/// This function is intended for when the result vector's length is known at
/// runtime, i.e., none of the inputs is a static vector.
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

    if (comptime !((meta.isStaticVector(X) or meta.isNumeric(X)) and (meta.isStaticVector(Y) or meta.isNumeric(Y)))) {
        if (x_len != y_len)
            return vector.Error.DimensionMismatch;
    }

    var result = switch (comptime meta.vectorType(R)) {
        .static => R.init,
        .dense => try R.init(allocator, x_len),
        .sparse => try R.init(allocator, x_len, apply2NNZ(if (comptime meta.isSparseVector(X)) x.idx[0..x.nnz] else &.{}, if (comptime meta.isSparseVector(Y)) y.idx[0..y.nnz] else &.{})),
        .numeric => unreachable,
    };

    vecops.apply2IntoUnchecked(
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
    );

    return result;
}

/// Applies a binary into operation elementwise between an output and two input
/// vectors, or between an output vector, an input vector and an input numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs, or sparse outputs when both inputs
/// are not vectors. Any other form of memory overlap might yield incorrect
/// results.
///
/// For three static vectors, or a static output vector, an input static vector
/// and an input numeric, the function cannot return an error unless opInto can.
///
/// For two input sparse vectors, or an input sparse vector and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the vectors has a non-zero element.
///
/// ## Signature
/// ```zig
/// vector.apply2Into(*O, x: X, y: Y, opInto: OpInto) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `opInto` (`comptime anytype`): An ito binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the vectors do not have the same
///   length.
pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        (!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.vector.apply2Into: o must be a mutable one-item pointer to a vector, at least one of x or y must be a vector, the other must be a vector or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (opInto != numeric.addInto and opInto != numeric.subInto and opInto != numeric.mulInto and opInto != numeric.divInto)
        @compileError("zsl.vector.apply2Into: opInto must be zsl.numeric.addInto, zsl.numeric.subInto, zsl.numeric.mulInto or zsl.numeric.divInto, got\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    O = meta.Child(O);

    const o_len = if (comptime meta.isStaticVector(O)) O.len else o.len;
    const x_len = if (comptime meta.isVector(X)) (if (comptime meta.isStaticVector(X)) X.len else x.len) else o_len;
    const y_len = if (comptime meta.isVector(Y)) (if (comptime meta.isStaticVector(Y)) Y.len else y.len) else o_len;

    if (comptime meta.isStaticVector(O) and
        (meta.isStaticVector(X) or meta.isNumeric(X)) and
        (meta.isStaticVector(Y) or meta.isNumeric(Y)))
    {
        if (comptime o_len != x_len or o_len != y_len)
            if (comptime meta.isVector(X) and meta.isVector(Y))
                @compileError("zsl.vector.apply2Into: static vectors o, x and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n")
            else if (comptime meta.isVector(X))
                @compileError("zsl.vector.apply2Into: static vectors o and x must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n")
            else
                @compileError("zsl.vector.apply2Into: static vectors o and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (o_len != x_len or o_len != y_len)
            return vector.Error.DimensionMismatch;
    }

    if (comptime meta.isSparseVector(O)) {
        if (comptime meta.isStaticVector(X) or meta.isDenseVector(X) or meta.isStaticVector(Y) or meta.isDenseVector(Y))
            @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

        const nnz = apply2NNZ(
            if (comptime meta.isSparseVector(X)) x.idx[0..x.nnz] else &.{},
            if (comptime meta.isSparseVector(Y)) y.idx[0..y.nnz] else &.{},
        );

        if (o._dlen < nnz or o._ilen < nnz)
            return vector.Error.InsufficientSpace;
    }

    return apply2IntoUnchecked(o, x, y, opInto);
}

/// Applies a binary into operation elementwise between an output and two input
/// vectors, or between an output vector, an input vector and an input numeric,
/// without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs, or sparse outputs when both inputs
/// are not vectors. Any other form of memory overlap might yield incorrect
/// results.
///
/// For two input sparse vectors, or an input sparse vector and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the vectors has a non-zero element.
///
/// ## Signature
/// ```zig
/// vector.apply2IntoUnchecked(*O, x: X, y: Y, opInto: OpInto) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `opInto` (`comptime anytype`): An ito binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        (!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.vector.apply2IntoUnchecked: o must be a mutable one-item pointer to a vector, at least one of x or y must be a vector, the other must be a vector or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (opInto != numeric.addInto and opInto != numeric.subInto and opInto != numeric.mulInto and opInto != numeric.divInto)
        @compileError("zsl.vector.apply2IntoUnchecked: opInto must be zsl.numeric.addInto, zsl.numeric.subInto, zsl.numeric.mulInto or zsl.numeric.divInto, got\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    O = meta.Child(O);

    switch (comptime meta.vectorType(O)) {
        .static => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecsta_vecsta_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecsta_vecsta_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecsta_vecsta_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecsta_vecsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecsta_vecden_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecsta_vecden_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecsta_vecden_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecsta_vecden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecsta_vecspa_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecsta_vecspa_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecsta_vecspa_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecsta_vecspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecsta_num_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecsta_num_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecsta_num_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .dense => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecden_vecsta_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecden_vecsta_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecden_vecsta_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecden_vecsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecden_vecden_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecden_vecden_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecden_vecden_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecden_vecden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecden_vecspa_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecden_vecspa_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecden_vecspa_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecden_vecspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/vecden_num_vecsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .dense => return @import("apply2/vecden_num_vecden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .sparse => return @import("apply2/vecden_num_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .sparse => switch (comptime meta.vectorType(X)) {
            .sparse => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2/vecspa_vecspa_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => return @import("apply2/vecspa_vecspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.vector.apply2IntoUnchecked: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2/vecspa_num_vecspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .numeric => unreachable,
                else => @compileError("zsl.vector.apply2IntoUnchecked: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            else => @compileError("zsl.vector.apply2IntoUnchecked: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
        },
        .numeric => unreachable,
    }
}

fn apply2NNZ(x_idx: []const usize, y_idx: []const usize) usize {
    var intersections: usize = 0;
    var px: usize = 0;
    var py: usize = 0;

    while (px < x_idx.len and py < y_idx.len) {
        const ix = x_idx[px];
        const iy = y_idx[py];

        if (ix == iy) {
            intersections += 1;
            px += 1;
            py += 1;
        } else if (ix < iy) {
            px += 1;
        } else {
            py += 1;
        }
    }

    return x_idx.len + y_idx.len - intersections;
}
