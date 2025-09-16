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

const utils = @import("utils.zig");

pub fn LDLT(T: type, order: Order) type {
    return struct {
        p: matrix.Permutation(T),
        l: matrix.Triangular(T, .lower, .unit, order),
        d: matrix.Tridiagonal(T), // strictly D is block diagonal with 1x1 and 2x2 blocks, but since no BDiagonal type exists, we use Tridiagonal here

        pub fn init(p: [*]u32, l: anytype, d: anytype, n: u32) LDLT(Numeric(Child(@TypeOf(l))), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = n,
                    .direction = .backward,
                    .flags = .{ .owns_data = true },
                },
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
                    .osize = n,
                    .offset = 0,
                    .sdoffset = (n - 1) + n,
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
        p: matrix.Permutation(T),
        u: matrix.Triangular(T, .upper, .unit, order),
        d: matrix.Tridiagonal(T), // strictly D is block diagonal with 1x1 and 2x2 blocks, but since no BDiagonal type exists, we use Tridiagonal here

        pub fn init(p: [*]u32, u: anytype, d: anytype, n: u32) UDUT(Numeric(Child(@TypeOf(u))), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = n,
                    .direction = .forward,
                    .flags = .{ .owns_data = true },
                },
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
                    .osize = n,
                    .offset = 0,
                    .sdoffset = (n - 1) + n,
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

pub fn ldlt(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !LDLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
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
        .symmetric => {
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

                var d: matrix.Tridiagonal(Numeric(A)) = try .init(allocator, n);
                errdefer d.deinit(allocator);

                var p: vector.Vector(u32) = try vector.Vector(u32).init(allocator, n);
                errdefer p.deinit(allocator);

                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    p.data[i] = i;
                }

                try lconvert(orderOf(A), n, l.data, d.data, p.data, ipiv.data, ctx);

                lpermute(n, p.data, ipiv.data);

                return .init(p.data, l.data, d.data, n);
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

                var d: matrix.Tridiagonal(Numeric(A)) = try .init(allocator, n);
                errdefer d.deinit(allocator);

                var p: vector.Vector(u32) = try vector.Vector(u32).init(allocator, n);
                errdefer p.deinit(allocator);

                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    p.data[i] = i;
                }

                try lconvert(orderOf(A), n, l.data, d.data, p.data, ipiv.data, ctx);

                lpermute(n, p.data, ipiv.data);

                return .init(p.data, l.data, d.data, n);
            }
        },
        .hermitian => @compileError("ldlt: hermitian matrices not implemented yet"),
        else => unreachable,
    }
}

/// Convert the output a of `sytrf` or `hetrf` into the l and d matrices of the
/// LDLT factorization. On output, `a` is `l` and `d` is overwritten. `p` is
/// used as workspace for the permutation.
fn lconvert(comptime order: Order, n: u32, a: anytype, d: anytype, p: [*]u32, ipiv: [*]i32, ctx: anytype) !void {
    var i: u32 = 0;
    while (i < n) : (i += 1) {
        if (ipiv[i] > 0) {
            // 1x1 pivot

            // Update p
            var tmp: u32 = p[i];
            p[i] = p[types.scast(u32, ipiv[i] - 1)];
            p[types.scast(u32, ipiv[i] - 1)] = tmp;

            // Extract d
            d[i + (n - 1)] = a[utils.cindex(order, i, i, n)]; // diagonal
            a[utils.cindex(order, i, i, n)] = undefined;

            if (i < n - 1) { // zero out under d11
                d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
            }

            // Permute rows of l
            while (p[i] != i) {
                const jrow: u32 = p[i];
                swap_rows_prefix(n, a, i, jrow, i);

                tmp = p[i];
                p[i] = p[jrow];
                p[jrow] = tmp;
            }
        } else {
            // 2x2 pivot

            // Update p
            var tmp: u32 = p[i + 1];
            p[i + 1] = p[types.scast(u32, -ipiv[i] - 1)];
            p[types.scast(u32, -ipiv[i] - 1)] = tmp;

            // Extract d
            d[i + (n - 1)] = a[utils.cindex(order, i, i, n)]; // diagonal (a11)
            d[i] = a[utils.cindex(order, i + 1, i, n)]; // subdiagonal (a21)
            d[i + (n - 1) + n] = try ops.copy(a[utils.cindex(order, i + 1, i, n)], ctx); // superdiagonal (a12)
            d[i + 1 + (n - 1)] = a[utils.cindex(order, i + 1, i + 1, n)]; // diagonal (a22)
            a[utils.cindex(order, i, i, n)] = undefined; // a11
            a[utils.cindex(order, i + 1, i, n)] = try constants.zero(Child(@TypeOf(a)), ctx); // a21
            a[utils.cindex(order, i + 1, i + 1, n)] = undefined; // a22

            i += 1;

            // Permute rows of l
            while (p[i] != i) {
                const jrow: u32 = p[i];
                swap_rows_prefix(n, a, i, jrow, i - 1);

                tmp = p[i];
                p[i] = p[jrow];
                p[jrow] = tmp;
            }

            if (i < n - 1) { // zero out under d22
                d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
            }
        }
    }
}

fn swap_rows_prefix(n: u32, a: anytype, r1: u32, r2: u32, upto_col: u32) void {
    var j: u32 = 0;
    while (j < upto_col) : (j += 1) {
        const tmp = a[r1 + j * n];
        a[r1 + j * n] = a[r2 + j * n];
        a[r2 + j * n] = tmp;
    }
}

fn lpermute(n: u32, p: [*]u32, ipiv: [*]i32) void {
    var i: u32 = 0;
    while (i < n) : (i += 1) {
        if (ipiv[i] > 0) {
            // 1x1 pivot
            const tmp: u32 = p[i];
            p[i] = p[types.scast(u32, ipiv[i] - 1)];
            p[types.scast(u32, ipiv[i] - 1)] = tmp;
        } else {
            // 2x2 pivot
            const tmp: u32 = p[i + 1];
            p[i + 1] = p[types.scast(u32, -ipiv[i] - 1)];
            p[types.scast(u32, -ipiv[i] - 1)] = tmp;

            i += 1;
        }
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
        .symmetric => {
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
        .hermitian => @compileError("udut: hermitian matrices not implemented yet"),
        else => unreachable,
    }
}
