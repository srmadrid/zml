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
        l: matrix.Triangular(T, .lower, .non_unit, order),

        pub fn init(l: anytype, n: u32) LLT(Numeric(Child(@TypeOf(l))), order) {
            return .{
                .l = .{
                    .data = l,
                    .rows = n,
                    .cols = n,
                    .ld = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *LLT(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self.l.data[0 .. self.l.rows * self.l.cols]);

            self.* = undefined;
        }
    };
}

pub fn UTU(T: type, order: Order) type {
    return struct {
        u: matrix.Triangular(T, .upper, .non_unit, order),

        pub fn init(u: anytype, n: u32) UTU(Numeric(Child(@TypeOf(u))), order) {
            return .{
                .u = .{
                    .data = u,
                    .rows = n,
                    .cols = n,
                    .ld = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *UTU(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self.u.data[0 .. self.u.rows * self.u.cols]);

            self.* = undefined;
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
                return error.SingularMatrix;

            return .init(l.data, n);
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

            return .init(u.data, n);
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
