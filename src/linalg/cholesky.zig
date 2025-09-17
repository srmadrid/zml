const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const orderOf = types.orderOf;
const Numeric = types.Numeric;
const Child = types.Child;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");

const linalg = @import("../linalg.zig");

pub fn LLT(T: type, order: Order) type {
    return struct {
        _size: u32,
        _l: [*]T,

        pub fn init(size: u32, a: [*]T) LLT(T, order) {
            return .{
                ._size = size,
                ._l = a,
            };
        }

        pub fn deinit(self: *LLT(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._l[0 .. self._size * self._size]);

            self.* = undefined;
        }

        pub fn l(self: *const LLT(T, order)) matrix.Triangular(T, .lower, .non_unit, order) {
            return .{
                .data = self._l,
                .rows = self._size,
                .cols = self._size,
                .ld = self._size,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn reconstruct(self: *const LLT(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = L * L^T
            @compileError("LLT.reconstruct not implemented yet");
        }
    };
}

pub fn UTU(T: type, order: Order) type {
    return struct {
        _size: u32,
        _u: [*]T,

        pub fn init(n: u32, a: [*]T) UTU(T, order) {
            return .{
                ._size = n,
                ._u = a,
            };
        }

        pub fn deinit(self: *UTU(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._u[0 .. self._size * self._size]);

            self.* = undefined;
        }

        pub fn u(self: *const UTU(T, order)) matrix.Triangular(T, .upper, .non_unit, order) {
            return .{
                .data = self._u,
                .rows = self._size,
                .cols = self._size,
                .ld = self._size,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn reconstruct(self: *const UTU(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = U * U^T
            @compileError("UTU.reconstruct not implemented yet");
        }
    };
}

pub fn llt(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !LLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("llt: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.llt not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .symmetric, .hermitian => {
            const n: u32 = a.size;

            var l = if (comptime types.uploOf(A) == .lower) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
            errdefer l.deinit(allocator);

            const info: i32 = try linalg.lapack.potrf(
                types.orderOf(A),
                .lower,
                types.scast(i32, n),
                l.data,
                types.scast(i32, l.ld),
                ctx,
            );

            if (info != 0)
                return linalg.Error.SingularMatrix;

            return .init(n, l.data);
        },
        else => unreachable,
    }
}

pub fn utu(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !UTU(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("utu: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.utu not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .symmetric, .hermitian => {
            const n: u32 = a.size;

            var u = if (comptime types.uploOf(A) == .upper) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
            errdefer u.deinit(allocator);

            const info: i32 = try linalg.lapack.potrf(
                types.orderOf(A),
                .upper,
                types.scast(i32, n),
                u.data,
                types.scast(i32, u.ld),
                ctx,
            );

            if (info != 0)
                return error.SingularMatrix;

            return .init(n, u.data);
        },
        else => unreachable,
    }
}

pub fn cholesky(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !(if (types.uploOf(@TypeOf(a)) == .lower) LLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) else UTU(Numeric(@TypeOf(a)), orderOf(@TypeOf(a)))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("cholesky: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.cholesky not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.uploOf(A)) {
        .lower => return llt(allocator, a, ctx),
        .upper => return utu(allocator, a, ctx),
    }
}
