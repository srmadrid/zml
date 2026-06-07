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

/// Applies a binary into operation elementwise between an output and two input
/// vectors, or between an output vector, an input vector and an input numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
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
        @compileError("zsl.vector.apply2Into: o must be a mutable one-itme pointer to a vector, at least one of x or y must be a vector, the other must be a vector or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

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
                @compileError("zsl.vector.apply2: static vectors o, x and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n")
            else if (comptime meta.isVector(X))
                @compileError("zsl.vector.apply2: static vectors o and x must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n")
            else
                @compileError("zsl.vector.apply2: static vectors o and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (o_len != x_len or o_len != y_len)
            return vector.Error.DimensionMismatch;
    }

    switch (comptime meta.vectorType(O)) {
        .static => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/ststst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/ststde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/ststsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/ststnu.zig").apply2Into(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/stdest.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/stdede.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/stdesp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/stdenu.zig").apply2Into(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/stspst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/stspde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/stspsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/stspnu.zig").apply2Into(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/stnust.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/stnude.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/stnusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .dense => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/destst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/destde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/destsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/destnu.zig").apply2Into(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/dedest.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/dedede.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/dedesp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/dedenu.zig").apply2Into(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/despst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/despde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/despsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/despnu.zig").apply2Into(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2/denust.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2/denude.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2/denusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .sparse => switch (comptime meta.vectorType(X)) {
            .sparse => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2/spspsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2/spspnu.zig").apply2Into(o, x, y, opInto),
                else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2/spnusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
                else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
        },
        .numeric => unreachable,
    }
}
