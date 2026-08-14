const std = @import("std");

const meta = @import("../meta.zig");

const vector = @import("../vector.zig");

const linalg = @import("../linalg.zig");

pub fn Normalize(comptime X: type, comptime norm_order: linalg.NormOrder(meta.Numeric(X))) type {
    comptime if (!meta.isVector(X))
        @compileError("zsl.linalg.Normalize: X must be a vector type, got\n\tX = " ++
            @typeName(X) ++ "\n");

    return vector.Div(X, linalg.Norm(X, norm_order));
}

/// Normalizes a static vector, `x / ‖x‖`, such that `‖x‖ = 1` according to a
/// chosen norm.
///
/// ## Signature
/// ```zig
/// linalg.normalize(x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) linalg.Normalize(X, order)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The vector.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `linalg.Normalize(@TypeOf(x), order)`: The normalized vector, `x / ‖x‖`.
pub fn normalize(x: anytype, comptime norm_order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Normalize(@TypeOf(x), norm_order) {
    const X: type = @TypeOf(x);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X))
        @compileError("zsl.linalg.normalize: the result cannot be a heap-allocated vector type, i.e., x must be a static vector, got\n\tx: " ++
            @typeName(X) ++ "\n\tresult: " ++ @typeName(linalg.Normalize(X, norm_order)) ++ "\nFor these inputs use zsl.linalg.normalizeAlloc instead.");

    return vector.div(x, linalg.norm(x, norm_order));
}

/// Normalizes a vector, `x / ‖x‖`, such that `‖x‖ = 1` according to a chosen
/// norm, dynamically allocating memory for the result.
///
/// ## Signature
/// ```zig
/// linalg.normalizeAlloc(allocator: std.mem.Allocator, x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) linalg.Normalize(X, order)
/// ```
///
/// ## Arguments
/// * allocator (`std.mem.Allocator): The allocator to use for memory`
///   allocations.
/// * `x` (`anytype`): The vector.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `linalg.Normalize(@TypeOf(x), order)`: The normalized vector, `x / ‖x‖`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn normalizeAlloc(allocator: std.mem.Allocator, x: anytype, comptime norm_order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) !linalg.Normalize(@TypeOf(x), norm_order) {
    return vector.divAlloc(allocator, x, linalg.norm(x, norm_order));
}

/// Normalizes a vector, `x / ‖x‖`, such that `‖x‖ = 1` according to a chosen
/// norm, storing the result in an output vector.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted.
///
/// ## Signature
/// ```zig
/// linalg.normalizeInto(o: *O, x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector.
/// * `x` (`anytype`): The input vector.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.Error.DimensionMismatch`: If the vectors do not have the same
///   length.
pub fn normalizeInto(o: anytype, x: anytype, comptime norm_order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or !meta.isVector(X))
        @compileError("zsl.linalg.normalizeInto: o must be a mutable one-item pointer to a vector, and x must be a vector, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    const o_len = if (comptime meta.isStaticVector(O)) O.len else o.len;
    const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;

    if (comptime meta.isStaticVector(O) and meta.isStaticVector(X)) {
        if (comptime o_len != x_len)
            @compileError("zsl.linalg.normalizeInto: static vectors o and x must have equal compile time lengths, got\n\to: *" ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");
    } else {
        if (o_len != x_len)
            return linalg.Error.DimensionMismatch;
    }

    return vector.divIntoUnchecked(o, x, linalg.norm(x, norm_order));
}

/// Normalizes a vector, `x / ‖x‖`, such that `‖x‖ = 1` according to a chosen
/// norm, storing the result in an output vector, and without performing any
/// length checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted.
///
/// ## Signature
/// ```zig
/// linalg.normalizeInto(o: *O, x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector.
/// * `x` (`anytype`): The input vector.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `void`
pub fn normalizeIntoUnchecked(o: anytype, x: anytype, comptime norm_order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or !meta.isVector(X))
        @compileError("zsl.linalg.normalizeIntoUnchecked: o must be a mutable one-item pointer to a vector, and x must be a vector, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    return vector.divIntoUnchecked(o, x, linalg.norm(x, norm_order));
}
