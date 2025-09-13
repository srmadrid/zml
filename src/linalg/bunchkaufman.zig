const std = @import("std");
const options = @import("options");

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

pub fn LDLT(T: type, order: Order) type {
    return struct {
        l: matrix.Triangular(T, .lower, .unit, order),
        d: matrix.Tridiagonal(T), // strictly D is block diagonal with 1x1 and 2x2 blocks, but since no BDiagonal type exists, we use Tridiagonal here
        p: matrix.Permutation(T),

        pub fn init(l: anytype, d: anytype, p: [*]u32, n: u32) LDLT(Numeric(Child(@TypeOf(l))), order) {
            return .{
                .l = .{
                    .data = l,
                    .rows = n,
                    .cols = n,
                    .ld = n,
                    .flags = .{ .owns_data = true },
                },
                .d = .{
                    .data = d,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
                .p = .{
                    .data = p,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *LDLT(T, order), allocator: std.mem.Allocator) void {
            self.l.deinit(allocator);
            self.d.deinit(allocator);
            self.p.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn UDUT(T: type, order: Order) type {
    return struct {
        u: matrix.Triangular(T, .upper, .unit, order),
        d: matrix.Tridiagonal(T), // strictly D is block diagonal with 1x1 and 2x2 blocks, but since no BDiagonal type exists, we use Tridiagonal here
        p: matrix.Permutation(T),

        pub fn init(u: anytype, d: anytype, p: [*]u32, n: u32) UDUT(Numeric(Child(@TypeOf(u))), order) {
            return .{
                .u = .{
                    .data = u,
                    .rows = n,
                    .cols = n,
                    .ld = n,
                    .flags = .{ .owns_data = true },
                },
                .d = .{
                    .data = d,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
                .p = .{
                    .data = p,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *UDUT(T, order), allocator: std.mem.Allocator) void {
            self.u.deinit(allocator);
            self.d.deinit(allocator);
            self.p.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn ldlt(allocator: std.mem.Allocator, a: anytype, ctx: anytype) ![*]Numeric(@TypeOf(a)) { // !LDLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("ldlt: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.ldlt not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .symmetric, .hermitian => {
            const n: u32 = a.size;

            var info: i32 = 0;
            if (comptime options.link_lapacke != null) {
                var l = if (comptime types.uploOf(A) == .lower) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer l.deinit(allocator);

                var ipiv: vector.Vector(i32) = try vector.Vector(i32).init(allocator, n);
                defer ipiv.deinit(allocator);

                info = try linalg.lapack.sytrf( // LAPACKE version does not need work
                    types.orderOf(A),
                    .lower,
                    types.scast(i32, n),
                    l.data,
                    types.scast(i32, l.ld),
                    ipiv.data,
                    undefined,
                    undefined,
                    ctx,
                );

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return l.data;
            } else {
                var lwork: i32 = -1;
                _ = try linalg.lapack.sytrf( // work size query
                    types.orderOf(A),
                    .lower,
                    types.scast(i32, n),
                    @as([*]f64, undefined),
                    types.scast(i32, a.ld),
                    undefined,
                    @as([*]i32, @ptrCast(&lwork)),
                    -1,
                    ctx,
                );

                var l = if (comptime types.uploOf(A) == .lower) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer l.deinit(allocator);

                var ipiv: vector.Vector(i32) = try vector.Vector(i32).init(allocator, n);
                defer ipiv.deinit(allocator);

                var work: vector.Vector(Numeric(A)) = try vector.Vector(Numeric(A)).init(allocator, types.scast(u32, lwork));
                defer work.deinit(allocator);

                info = try linalg.lapack.sytrf(
                    types.orderOf(A),
                    .lower,
                    types.scast(i32, n),
                    l.data,
                    types.scast(i32, l.ld),
                    ipiv.data,
                    work.data,
                    types.scast(i32, work.len),
                    ctx,
                );

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return l.data;
            }
        },
        else => unreachable,
    }
}

pub fn udut(allocator: std.mem.Allocator, a: anytype, ctx: anytype) ![*]Numeric(@TypeOf(a)) { // !UDUT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("udut: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.udut not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .symmetric, .hermitian => {
            const n: u32 = a.size;

            var info: i32 = 0;
            if (comptime options.link_lapacke != null) {
                var u = if (comptime types.uploOf(A) == .upper) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer u.deinit(allocator);

                var ipiv: vector.Vector(i32) = try vector.Vector(i32).init(allocator, n);
                defer ipiv.deinit(allocator);

                info = try linalg.lapack.sytrf( // LAPACKE version does not need work
                    types.orderOf(A),
                    .upper,
                    types.scast(i32, n),
                    u.data,
                    types.scast(i32, u.ld),
                    ipiv.data,
                    undefined,
                    undefined,
                    ctx,
                );

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return u.data;
            } else {
                var lwork: i32 = -1;
                _ = try linalg.lapack.sytrf( // work size query
                    types.orderOf(A),
                    .upper,
                    types.scast(i32, n),
                    @as([*]f64, undefined),
                    types.scast(i32, a.ld),
                    undefined,
                    @as([*]i32, @ptrCast(&lwork)),
                    -1,
                    ctx,
                );

                var u = if (comptime types.uploOf(A) == .upper) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer u.deinit(allocator);

                var ipiv: vector.Vector(i32) = try vector.Vector(i32).init(allocator, n);
                defer ipiv.deinit(allocator);

                var work: vector.Vector(Numeric(A)) = try vector.Vector(Numeric(A)).init(allocator, types.scast(u32, lwork));
                defer work.deinit(allocator);

                info = try linalg.lapack.sytrf(
                    types.orderOf(A),
                    .upper,
                    types.scast(i32, n),
                    u.data,
                    types.scast(i32, u.ld),
                    ipiv.data,
                    work.data,
                    types.scast(i32, work.len),
                    ctx,
                );

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return u.data;
            }
        },
        else => unreachable,
    }
}
