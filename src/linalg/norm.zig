const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const linalg = @import("../linalg.zig");

pub fn Norm(comptime X: type, comptime order: linalg.NormOrder(meta.Numeric(X))) type {
    comptime if (!meta.isVector(X) and !meta.isMatrix(X))
        @compileError("zsl.linalg.Norm: X must be a vector or matrix type, got\n\tX = " ++
            @typeName(X) ++ "\n");

    return switch (comptime order) {
        .l1 => numeric.Add(
            numeric.Abs(meta.Numeric(X)),
            numeric.Abs(meta.Numeric(X)),
        ),
        .l2, .frobenius => numeric.Sqrt(
            numeric.Add(
                numeric.Abs(meta.Numeric(X)),
                numeric.Abs(meta.Numeric(X)),
            ),
        ),
        .inf => numeric.Max(
            numeric.Abs(meta.Numeric(X)),
            numeric.Abs(meta.Numeric(X)),
        ),
        .p => numeric.Pow(
            numeric.Add(
                numeric.Pow(numeric.Abs(meta.Numeric(X)), meta.Numeric(X)),
                numeric.Pow(numeric.Abs(meta.Numeric(X)), meta.Numeric(X)),
            ),
            meta.Numeric(X),
        ),
    };
}

pub fn NormOrder(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.NormOrder: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return union(enum) {
        l1,
        l2,
        frobenius,
        inf,
        p: N,
    };
}

/// Computes the norm of a vector or a matrix, ‖x‖. The following norm orders
/// are supported for vectors:
/// * ‖x‖₁ = ∑ᵢ |xᵢ|
/// * ‖x‖₂ = √(∑ᵢ |xᵢ|²)
/// * ‖x‖_F = √(∑ᵢ |xᵢ|²)
/// * ‖x‖_∞ = maxᵢ |xᵢ|
/// * ‖x‖ₚ = (∑ᵢ |xᵢ|ᵖ)¹^/ᵖ
/// The following norm orders are supported for matrices:
/// * ‖x‖₁ = maxⱼ ∑ᵢ |xᵢⱼ|
/// * ‖x‖_F = √(∑ᵢ∑ⱼ |xᵢⱼ|²)
/// * ‖x‖_∞ = maxᵢ ∑ⱼ |xᵢⱼ|
///
/// ## Signature
/// ```zig
/// linalg.norm(x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) linalg.Norm(X, order)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The matrix.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `linalg.Norm(@TypeOf(x), order)`: The norm ‖x‖.
pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X: type = @TypeOf(x);

    if (comptime meta.isMatrix(X) and (meta.isDenseMatrix(X) or meta.isSparseMatrix(X)) and order == .l2)
        @compileError("zsl.linalg.norm: the l2-norm of heap-allocated matrices requires internal temporary allocations. For these inputs use zsl.linalg.normAlloc instead.");

    if (comptime meta.isMatrix(X) and order == .p)
        @compileError("zsl.linalg.norm: the p-norm of matrices is not supported.");

    switch (comptime meta.domain(X)) {
        .vector => switch (comptime meta.vectorType(X)) {
            .static => return @import("norm/vecsta.zig").norm(x, order),
            // .dense => return @import("norm/de.zig").norm(x, norm_type),
            // .sparse => return @import("norm/sp.zig").norm(x, norm_type),
            else => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n"),
            .numeric => unreachable,
        },
        .matrix => switch (comptime meta.matrixType(X)) {
            .general_dense => return @import("norm/matgenden.zig").norm(x, order),
            else => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n"), // Only dense and sparse matrices need allocator, static ones can use stack allocated work array
        },
        else => unreachable,
    }
}

/// Computes the norm of a vector or a matrix, ‖x‖, allocating when necessary or
/// beneficial. The following norm orders are supported for vectors:
/// * ‖x‖₁ = ∑ᵢ |xᵢ|
/// * ‖x‖₂ = √(∑ᵢ |xᵢ|²)
/// * ‖x‖_F = √(∑ᵢ |xᵢ|²)
/// * ‖x‖_∞ = maxᵢ |xᵢ|
/// * ‖x‖ₚ = (∑ᵢ |xᵢ|ᵖ)¹^/ᵖ
/// The following norm orders are supported for matrices:
/// * ‖x‖₁ = maxⱼ ∑ᵢ |xᵢⱼ|
/// * ‖x‖₂ = maxᵢ σᵢ(x) = maxᵢ √(λᵢ(Aᴴ A))
/// * ‖x‖_F = √(∑ᵢ∑ⱼ |xᵢⱼ|²)
/// * ‖x‖_∞ = maxᵢ ∑ⱼ |xᵢⱼ|
///
/// ## Signature
/// ```zig
/// linalg.normAlloc(allocator: std.mem.Allocator, x: X, comptime order: linalg.NormOrder(meta.Numeric(X))) !linalg.Norm(X, order)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The matrix.
/// * `order` (`comptime linalg.NormOrder(meta.Numeric(X))`): the norm order.
///
/// ## Returns
/// `linalg.Norm(@TypeOf(x), order)`: The norm ‖x‖.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn normAlloc(allocator: std.mem.Allocator, x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) !linalg.Norm(@TypeOf(x), order) {
    const X: type = @TypeOf(x);

    if (comptime meta.isMatrix(X) and order == .p)
        @compileError("zsl.linalg.normAlloc: the p-norm of matrices is not supported.");

    switch (comptime meta.domain(X)) {
        .vector => switch (comptime meta.vectorType(X)) {
            .static => return @import("norm/vecsta.zig").norm(x, order),
            // .dense => return @import("norm/de.zig").norm(x, norm_type),
            // .sparse => return @import("norm/sp.zig").norm(x, norm_type),
            else => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n"),
            .numeric => unreachable,
        },
        .matrix => switch (comptime meta.matrixType(X)) {
            .general_dense => return @import("norm/matgenden.zig").normAlloc(allocator, x, order),
            else => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n"), // Only dense and sparse matrices need allocator, static ones can use stack allocated work array
        },
        else => unreachable,
    }
}
