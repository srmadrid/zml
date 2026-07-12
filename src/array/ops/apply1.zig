const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const array = @import("../../array.zig");

const arrutils = @import("../utils.zig");

pub fn Apply1(comptime X: type, comptime op: anytype) type {
    const Op = @TypeOf(op);
    const opinfo = @typeInfo(Op);

    comptime if (!meta.isArray(X) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 1)
        @compileError("zsl.array.Apply1: X must be an array type, and op must be a function of one argument, got\n\tX = " ++
            @typeName(X) ++ "\n\tOp = " ++ @typeName(Op) ++ "\n");

    comptime var R = meta.ReturnTypeFromInputs(op, &.{meta.Numeric(X)});
    const rinfo = @typeInfo(R);
    if (rinfo == .error_union)
        R = rinfo.error_union.payload;

    comptime if (!meta.isNumeric(R))
        @compileError("zsl.array.Apply1: calling op with an argument of type " ++ @typeName(meta.Numeric(X)) ++ " must return a numeric, got\n\tR = " ++ @typeName(R) ++ "\n");

    switch (comptime meta.arrayType(X)) {
        .static => return array.Static(X.shape[0..X.ndim], R, X.storage_order),
        .dense => return array.Dense(R),
        .sparse => {
            const r = comptime op(numeric.zero(meta.Numeric(X)));

            return if (comptime numeric.eq(r, 0))
                array.Sparse(R, X.storage_order)
            else
                array.Dense(R);
        },
        else => unreachable,
    }
}

/// Applies a unary operation elementwise to an array.
///
/// This function is intended for when the result arrays's shape is known at
/// compile time, i.e., the input is a static array.
///
/// ## Signature
/// ```zig
/// array.apply1(x: X, op: Op) array.Apply1(X, op)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The operand.
/// * `op` (`comptime anytype`): A unary numeric function to apply elementwise
///   to `x`.
///
/// ## Returns
/// `array.Apply1(@TypeOf(x), op)`: The result of the operation.
pub fn apply1(x: anytype, comptime op: anytype) array.Apply1(@TypeOf(x), op) {
    const X: type = @TypeOf(x);
    const Op: type = @TypeOf(op);
    const R: type = array.Apply1(X, op);

    if (comptime meta.isDenseArray(R) or meta.isSparseArray(R))
        @compileError("zsl.array.apply1: the result cannot be a heap-allocated array type, i.e., the input must be a static array, got\n\tx: " ++
            @typeName(X) ++ "\n\top: " ++ @typeName(Op) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.array.apply1Alloc instead.");

    const NR = meta.ReturnTypeFromInputs(op, &.{meta.Numeric(X)});
    const nrinfo = @typeInfo(NR);

    const opInto = struct {
        pub fn opInto(o: anytype, nx: anytype) if (nrinfo == .error_union) anyerror!void else void {
            o.* = numeric.cast(
                meta.Numeric(meta.Numeric(R)),
                if (comptime nrinfo == .error_union)
                    try op(nx)
                else
                    op(nx),
            );
        }
    }.opInto;

    var result = R.init;

    array.apply1IntoUnchecked(
        &result,
        x,
        opInto,
    );

    return result;
}

/// Applies a unary operation elementwise to an array, dynamically allocating
/// memory for the result.
///
/// This function is intended for when the result array's shape is known at
/// runtime, i.e., the input is not a static array.
///
/// For a sparse array, the operation is only applied to the indices where it
/// has a non-zero element.
///
/// ## Signature
/// ```zig
/// array.apply1Alloc(allocator: std.mem.Allocator, x: X, op: Op) !array.Apply1(X, op)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The operand.
/// * `op` (`comptime anytype`): A unary numeric function to apply elementwise
///   to `x`.
///
/// ## Returns
/// `array.Apply1(@TypeOf(x), op)`: The result of the operation.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn apply1Alloc(allocator: std.mem.Allocator, x: anytype, comptime op: anytype) !array.Apply1(@TypeOf(x), op) {
    const X: type = @TypeOf(x);
    const R: type = array.Apply1(X, op);

    const NR = meta.ReturnTypeFromInputs(op, &.{meta.Numeric(X)});
    const nrinfo = @typeInfo(NR);

    const opInto = struct {
        pub fn opInto(o: anytype, nx: anytype) if (nrinfo == .error_union) anyerror!void else void {
            o.* = numeric.cast(
                meta.Numeric(meta.Numeric(R)),
                if (comptime nrinfo == .error_union)
                    try op(nx)
                else
                    op(nx),
            );
        }
    }.opInto;

    var result = switch (comptime meta.arrayType(R)) {
        .static => R.init,
        .dense => try R.init(
            allocator,
            x.shape[0..x.ndim],
            if (comptime meta.isDenseArray(X))
                arrutils.orderFromStrides(x.strides[0..x.ndim])
            else
                X.storage_order,
        ),
        .sparse => try R.initLike(allocator, x),
        else => unreachable,
    };

    array.apply1IntoUnchecked(
        &result,
        x,
        opInto,
    );

    return result;
}

/// Applies a unary into operation elementwise between an output and an input
/// array.
///
/// Exact aliasing (in-place modification) between the output and the input is
/// permitted.
///
/// For two static arrays, the function cannot return an error unless opInto
/// can.
///
/// When the output and the input ae sparse arrays, the operation is only
/// applied to the indices where the input has a non-zero element.
///
/// ## Signature
/// ```zig
/// array.apply1Into(*O, x: X, opInto: OpInto) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The input operand.
/// * `opInto` (`comptime anytype`): An into unary numeric function to apply
///   elementwise to `o` and `x`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `array.Error.DimensionMismatch`: If the arrays do not have compatible
///   shapes.
pub fn apply1Into(o: anytype, x: anytype, comptime opInto: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isArray(meta.Child(O)) or
        !meta.isArray(X) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 2)
        @compileError("zsl.array.apply1Into: o must be a mutable one-item pointer to an array, x must be an array, and opInto must be a function of two arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (meta.isBuilderSparseArray(X))
        @compileError("zsl.array.apply1Into: builder array types are not allowed as inputs, got\n\tx: " ++
            @typeName(X) ++ "\n");

    O = meta.Child(O);

    const o_shape = if (comptime meta.isStaticArray(O)) O.shape[0..O.ndim] else o.shape[0..o.ndim];
    const x_shape = if (comptime meta.isStaticArray(X)) X.shape[0..X.ndim] else x.shape[0..x.ndim];

    const broadcast_shape = if (comptime meta.isStaticArray(O) and meta.isStaticArray(X))
        comptime arrutils.broadcastShapes(&.{ o_shape, x_shape }) catch
            @compileError("zsl.array.apply1Into: static arrays o and x must have compatible compile time shapes, got\n\to: *" ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n")
    else
        try arrutils.broadcastShapes(&.{ o_shape, x_shape });
    if (comptime meta.isStaticArray(O) and meta.isStaticArray(X)) {
        if (comptime !std.mem.eql(usize, o_shape, broadcast_shape.shape[0..broadcast_shape.ndim]))
            @compileError("zsl.array.apply1Into: static arrays o and x must have compatible compile time shapes, got\n\to: *" ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");
    } else {
        if (!std.mem.eql(usize, o_shape, broadcast_shape.shape[0..broadcast_shape.ndim]))
            return array.Error.DimensionMismatch;
    }

    if (comptime meta.isSparseArray(O)) {
        if (comptime meta.isStaticArray(X) or meta.isDenseArray(X))
            @compileError("zsl.array.apply1Into: o cannot point to a sparse array if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

        if (comptime meta.isBuilderSparseArray(O)) {
            if (o._dlen < x.nnz or o._ilen < x.nnz)
                return array.Error.InsufficientSpace;
        } else {
            var cum_broadcast: usize = 1;
            var b_factors: [array.max_dimensions]usize = undefined;
            for (0..x.ndim) |level| {
                const dim = switch (X.storage_order) {
                    .f => x.ndim - 1 - level,
                    .c => level,
                };

                const x_dim_len = x.shape[dim];
                const o_dim_len = broadcast_shape.shape[dim];

                const factor = if (x_dim_len == 1 and o_dim_len > 1) o_dim_len else 1;
                cum_broadcast *= factor;
                b_factors[level] = cum_broadcast;
            }

            if (o._dlen < x.nnz * cum_broadcast)
                return array.Error.InsufficientSpace;

            var x_nodes_at_level: usize = 1;
            for (0..x.ndim) |level| {
                x_nodes_at_level = if (level == x.ndim - 1)
                    x.nnz
                else
                    x.ptr[level][x_nodes_at_level];

                const required_ilen = x_nodes_at_level * b_factors[level];
                if (o._ilen[level] < required_ilen)
                    return array.Error.InsufficientSpace;

                if (level > 0) {
                    const required_plen = (x_nodes_at_level * b_factors[level - 1]) + 1;
                    if (o._plen[level] < required_plen)
                        return array.Error.InsufficientSpace;
                }
            }
        }
    }

    return apply1IntoUnchecked(o, x, opInto);
}

/// Applies a unary into operation elementwise between an output and and input
/// array, without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted.
///
/// When the output and the input are sparse arrays, the operation is only
/// applied to the indices where the input has a non-zero element.
///
/// ## Signature
/// ```zig
/// array.apply1IntoUnchecked(*O, x: X, opInto: OpInto) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The input operand.
/// * `opInto` (`comptime anytype`): An into unary numeric function to apply
///   elementwise to `o` and `x`.
///
/// ## Returns
/// `void`
pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isArray(meta.Child(O)) or
        !meta.isArray(X) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 2)
        @compileError("zsl.array.apply1IntoUnchecked: o must be a mutable one-item pointer to an array, x must be an array, and opInto must be a function of two arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (meta.isBuilderSparseArray(X))
        @compileError("zsl.array.apply1IntoUnchecked: builder array types are not allowed as inputs, got\n\tx: " ++
            @typeName(X) ++ "\n");

    O = meta.Child(O);

    switch (comptime meta.arrayType(O)) {
        .static => switch (comptime meta.arrayType(X)) {
            .static => return @import("apply1/arrsta_arrsta.zig").apply1IntoUnchecked(o, x, opInto),
            .dense => return @import("apply1/arrsta_arrden.zig").apply1IntoUnchecked(o, x, opInto),
            .sparse => return @import("apply1/arrsta_arrspa.zig").apply1IntoUnchecked(o, x, opInto),
            else => unreachable,
        },
        .dense => switch (comptime meta.arrayType(X)) {
            .static => return @import("apply1/arrden_arrsta.zig").apply1IntoUnchecked(o, x, opInto),
            .dense => return @import("apply1/arrden_arrden.zig").apply1IntoUnchecked(o, x, opInto),
            .sparse => return @import("apply1/arrden_arrspa.zig").apply1IntoUnchecked(o, x, opInto),
            else => unreachable,
        },
        .sparse => switch (comptime meta.arrayType(X)) {
            .sparse => {
                comptime var no: meta.Numeric(O) = undefined;
                comptime opInto(&no, numeric.zero(meta.Numeric(X)));
                comptime if (numeric.ne(no, 0))
                    @compileError("zsl.array.apply1IntoUnchecked: opInto(&o, 0) must set o to zero when o and x are sparse arrays, got\n\to: *" ++
                        @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

                return @import("apply1/arrspa_arrspa.zig").apply1IntoUnchecked(o, x, opInto);
            },
            else => @compileError("zsl.array.apply1IntoUnchecked: o cannot point to a sparse array if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
        },
        .builder_sparse => switch (comptime meta.arrayType(X)) {
            .sparse => {
                comptime var no: meta.Numeric(O) = undefined;
                comptime opInto(&no, numeric.zero(meta.Numeric(X)));
                comptime if (numeric.ne(no, 0))
                    @compileError("zsl.array.apply1IntoUnchecked: opInto(&o, 0) must set o to zero when o and x are sparse arrays, got\n\to: *" ++
                        @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

                return @import("apply1/arrbuispa_arrspa.zig").apply1IntoUnchecked(o, x, opInto);
            },
            else => @compileError("zsl.array.apply1IntoUnchecked: o cannot point to a sparse array if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
        },
        .numeric => unreachable,
    }
}
