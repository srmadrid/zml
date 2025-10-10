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
        _size: u32,
        _ipiv: [*]i32, // 1-based indexing
        _ld: [*]T,
        _hermitian: bool,

        pub fn init(size: u32, a: [*]T, ipiv: [*]i32, hermitian: bool) LDLT(T, order) {
            return .{
                ._size = size,
                ._ipiv = ipiv,
                ._ld = a,
                ._hermitian = hermitian,
            };
        }

        pub fn deinit(self: *LDLT(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._ipiv[0..self._size]);
            allocator.free(self._ld[0 .. self._size * self._size]);

            self.* = undefined;
        }

        pub fn p(self: *const LDLT(T, order), allocator: std.mem.Allocator) !matrix.Permutation(T) {
            var _p: matrix.Permutation(T) = try .init(allocator, self._size);
            errdefer _p.deinit(allocator);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                _p.data[i] = i;
            }

            lpermute(self._size, _p.data, self._ipiv);

            return _p.transpose();
        }

        pub fn l(self: *const LDLT(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.Triangular(T, .lower, .unit, order) {
            var _l: matrix.Triangular(T, .lower, .unit, order) = try .init(allocator, self._size, self._size);
            errdefer _l.deinit(allocator);

            var pwork: [*]u32 = (try allocator.alloc(u32, self._size)).ptr;
            defer allocator.free(pwork[0..self._size]);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                pwork[i] = i;
            }

            try lconvert(order, self._size, self._ld, _l.data, null, pwork, self._ipiv, self._hermitian, ctx);

            return _l;
        }

        pub fn d(self: *const LDLT(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.Tridiagonal(T) {
            var _d: matrix.Tridiagonal(T) = try .init(allocator, self._size);
            errdefer _d.deinit(allocator);

            var pwork: [*]u32 = (try allocator.alloc(u32, self._size)).ptr;
            defer allocator.free(pwork[0..self._size]);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                pwork[i] = i;
            }

            try lconvert(order, self._size, self._ld, null, _d.data, pwork, self._ipiv, self._hermitian, ctx);

            return _d;
        }

        pub fn reconstruct(self: *const LDLT(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = P * L * D * L^T * P^T or A = P * L * D * L^H * P^T
            @compileError("LDLT.reconstruct not implemented yet");
        }
    };
}

pub fn UDUT(T: type, order: Order) type {
    return struct {
        _size: u32,
        _ipiv: [*]i32, // 1-based indexing
        _ud: [*]T,
        _hermitian: bool,

        pub fn init(size: u32, a: [*]T, ipiv: [*]i32, hermitian: bool) UDUT(T, order) {
            return .{
                ._size = size,
                ._ipiv = ipiv,
                ._ud = a,
                ._hermitian = hermitian,
            };
        }

        pub fn deinit(self: *UDUT(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._ipiv[0..self._size]);
            allocator.free(self._ud[0 .. self._size * self._size]);

            self.* = undefined;
        }

        pub fn p(self: *const UDUT(T, order), allocator: std.mem.Allocator) !matrix.Permutation(T) {
            var _p: matrix.Permutation(T) = try .init(allocator, self._size);
            errdefer _p.deinit(allocator);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                _p.data[i] = i;
            }

            upermute(self._size, _p.data, self._ipiv);

            return _p.transpose();
        }

        pub fn u(self: *const UDUT(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.Triangular(T, .upper, .unit, order) {
            var _u: matrix.Triangular(T, .upper, .unit, order) = try .init(allocator, self._size, self._size);
            errdefer _u.deinit(allocator);

            var pwork: [*]u32 = (try allocator.alloc(u32, self._size)).ptr;
            defer allocator.free(pwork[0..self._size]);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                pwork[i] = i;
            }

            try uconvert(order, self._size, self._ud, _u.data, null, pwork, self._ipiv, self._hermitian, ctx);

            return _u;
        }

        pub fn d(self: *const UDUT(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.Tridiagonal(T) {
            var _d: matrix.Tridiagonal(T) = try .init(allocator, self._size);
            errdefer _d.deinit(allocator);

            var pwork: [*]u32 = (try allocator.alloc(u32, self._size)).ptr;
            defer allocator.free(pwork[0..self._size]);

            var i: u32 = 0;
            while (i < self._size) : (i += 1) {
                pwork[i] = i;
            }

            try uconvert(order, self._size, self._ud, null, _d.data, pwork, self._ipiv, self._hermitian, ctx);

            return _d;
        }

        pub fn reconstruct(self: *const UDUT(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = P * U * D * U^T * P^T or A = P * U * D * U^H * P^T
            @compileError("UDUT.reconstruct not implemented yet");
        }
    };
}

pub fn ldlt(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !LDLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("ldlt: a must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.ldlt not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .symmetric, .hermitian => {
            const n: u32 = a.size;

            if (comptime options.link_lapacke != null) {
                var ld = if (comptime types.uploOf(A) == .lower) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer ld.deinit(allocator);

                const ipiv: []i32 = try allocator.alloc(i32, n);
                errdefer allocator.free(ipiv);

                var info: i32 = 0;
                if (comptime !types.isHermitianMatrix(A)) {
                    info = try linalg.lapack.sytrf( // LAPACKE version does not need work
                        types.orderOf(A),
                        .lower,
                        types.scast(i32, n),
                        ld.data,
                        types.scast(i32, ld.ld),
                        ipiv.ptr,
                        undefined,
                        undefined,
                        ctx,
                    );
                } else {
                    info = try linalg.lapack.hetrf( // LAPACKE version does not need work
                        types.orderOf(A),
                        .lower,
                        types.scast(i32, n),
                        ld.data,
                        types.scast(i32, ld.ld),
                        ipiv.ptr,
                        undefined,
                        undefined,
                        ctx,
                    );
                }

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return .init(n, ld.data, ipiv.ptr, types.isHermitianMatrix(A));
            } else {
                var lwork: i32 = -1;
                if (comptime !types.isHermitianMatrix(A)) {
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
                } else {
                    _ = try linalg.lapack.hetrf( // work size query
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
                }

                var ld = if (comptime types.uploOf(A) == .lower) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer ld.deinit(allocator);

                const ipiv: []i32 = try allocator.alloc(i32, n);
                errdefer allocator.free(ipiv);

                const work: []Numeric(A) = try allocator.alloc(Numeric(A), types.scast(u32, lwork));
                defer allocator.free(work);

                var info: i32 = 0;
                if (comptime !types.isHermitianMatrix(A)) {
                    info = try linalg.lapack.sytrf(
                        types.orderOf(A),
                        .lower,
                        types.scast(i32, n),
                        ld.data,
                        types.scast(i32, ld.ld),
                        ipiv.ptr,
                        work.ptr,
                        lwork,
                        ctx,
                    );
                } else {
                    info = try linalg.lapack.hetrf(
                        types.orderOf(A),
                        .lower,
                        types.scast(i32, n),
                        ld.data,
                        types.scast(i32, ld.ld),
                        ipiv.ptr,
                        work.ptr,
                        lwork,
                        ctx,
                    );
                }

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return .init(n, ld.data, ipiv.ptr, types.isHermitianMatrix(A));
            }
        },
        else => unreachable,
    }
}

/// Convert the output a of `sytrf` or `hetrf` into the l and d matrices of the
/// LDLT factorization. On output, `a` is `l` and `d` is overwritten. `p` is
/// used as workspace for the permutation.
fn lconvert(
    comptime order: Order,
    n: u32,
    a: anytype,
    l: anytype,
    d: anytype,
    p: [*]u32,
    ipiv: [*]i32,
    hermitian: bool,
    ctx: anytype,
) !void {
    var i: u32 = 0;
    while (i < n) : (i += 1) {
        if (ipiv[i] > 0) {
            // 1x1 pivot

            // Update p
            var tmp: u32 = p[i];
            p[i] = p[types.scast(u32, ipiv[i] - 1)];
            p[types.scast(u32, ipiv[i] - 1)] = tmp;

            // Extract d
            if (comptime @TypeOf(d) != @TypeOf(null)) {
                d[i + (n - 1)] = try ops.copy(a[utils.cindex(order, i, i, n)], ctx); // diagonal

                if (i < n - 1) { // zero out under d11
                    d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                    d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
                }
            }

            // Extract l
            if (comptime @TypeOf(l) != @TypeOf(null)) {
                var j: u32 = i + 1;
                while (j < n) : (j += 1) {
                    l[utils.cindex(order, j, i, n)] = try ops.copy(a[utils.cindex(order, j, i, n)], ctx);
                }
            }

            // Permute rows of l
            while (p[i] != i) {
                const jrow: u32 = p[i];

                if (comptime @TypeOf(l) != @TypeOf(null))
                    swap_rows_prefix(order, n, l, i, jrow, i);

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
            if (comptime @TypeOf(d) != @TypeOf(null)) {
                d[i + (n - 1)] = try ops.copy(a[utils.cindex(order, i, i, n)], ctx); // diagonal (a11)
                d[i] = try ops.copy(a[utils.cindex(order, i + 1, i, n)], ctx); // subdiagonal (a21)
                d[i + (n - 1) + n] = if (hermitian)
                    try ops.conj(a[utils.cindex(order, i + 1, i, n)], ctx)
                else
                    try ops.copy(a[utils.cindex(order, i + 1, i, n)], ctx); // superdiagonal (a12)
                d[i + 1 + (n - 1)] = try ops.copy(a[utils.cindex(order, i + 1, i + 1, n)], ctx); // diagonal (a22)
            }

            // Extract l
            if (comptime @TypeOf(l) != @TypeOf(null)) {
                l[utils.cindex(order, i + 1, i, n)] = try constants.zero(Child(@TypeOf(a)), ctx);

                var j: u32 = i + 2;
                while (j < n) : (j += 1) {
                    l[utils.cindex(order, j, i, n)] = try ops.copy(a[utils.cindex(order, j, i, n)], ctx);
                    l[utils.cindex(order, j, i + 1, n)] = try ops.copy(a[utils.cindex(order, j, i + 1, n)], ctx);
                }
            }

            i += 1;

            // Permute rows of l
            while (p[i] != i) {
                const jrow: u32 = p[i];

                if (comptime @TypeOf(l) != @TypeOf(null))
                    swap_rows_prefix(order, n, l, i, jrow, i - 1);

                tmp = p[i];
                p[i] = p[jrow];
                p[jrow] = tmp;
            }

            if (comptime @TypeOf(d) != @TypeOf(null)) {
                if (i < n - 1) { // zero out under d22
                    d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                    d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
                }
            }
        }
    }
}

// swap rows r1 <-> r2 across columns [0 .. upto_col-1]
fn swap_rows_prefix(comptime order: Order, n: u32, a: anytype, r1: u32, r2: u32, upto_col: u32) void {
    var j: u32 = 0;
    while (j < upto_col) : (j += 1) {
        const tmp = a[utils.cindex(order, r1, j, n)];
        a[utils.cindex(order, r1, j, n)] = a[utils.cindex(order, r2, j, n)];
        a[utils.cindex(order, r2, j, n)] = tmp;
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

pub fn udut(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !UDUT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("udut: a must be a symmetric or hermitian matrix, got " ++ @typeName(A));

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
                var ud = if (comptime types.uploOf(A) == .upper) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer ud.deinit(allocator);

                const ipiv: []i32 = try allocator.alloc(i32, n);
                errdefer allocator.free(ipiv);

                if (comptime !types.isHermitianMatrix(A)) {
                    info = try linalg.lapack.sytrf( // LAPACKE version does not need work
                        types.orderOf(A),
                        .upper,
                        types.scast(i32, n),
                        ud.data,
                        types.scast(i32, ud.ld),
                        ipiv.ptr,
                        undefined,
                        undefined,
                        ctx,
                    );
                } else {
                    info = try linalg.lapack.hetrf( // LAPACKE version does not need work
                        types.orderOf(A),
                        .upper,
                        types.scast(i32, n),
                        ud.data,
                        types.scast(i32, ud.ld),
                        ipiv.ptr,
                        undefined,
                        undefined,
                        ctx,
                    );
                }

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return .init(n, ud.data, ipiv.ptr, types.isHermitianMatrix(A));
            } else {
                var lwork: i32 = -1;
                if (comptime !types.isHermitianMatrix(A)) {
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
                } else {
                    _ = try linalg.lapack.hetrf( // work size query
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
                }

                var ud = if (comptime types.uploOf(A) == .upper) try a.copy(allocator, ctx) else try a.copyInverseUplo(allocator, ctx);
                errdefer ud.deinit(allocator);

                const ipiv: []i32 = try allocator.alloc(i32, n);
                errdefer allocator.free(ipiv);

                const work: []Numeric(A) = try allocator.alloc(Numeric(A), types.scast(u32, lwork));
                defer allocator.free(work);

                if (comptime !types.isHermitianMatrix(A)) {
                    info = try linalg.lapack.sytrf(
                        types.orderOf(A),
                        .upper,
                        types.scast(i32, n),
                        ud.data,
                        types.scast(i32, ud.ld),
                        ipiv.ptr,
                        work.ptr,
                        lwork,
                        ctx,
                    );
                } else {
                    info = try linalg.lapack.hetrf(
                        types.orderOf(A),
                        .upper,
                        types.scast(i32, n),
                        ud.data,
                        types.scast(i32, ud.ld),
                        ipiv.ptr,
                        work.ptr,
                        lwork,
                        ctx,
                    );
                }

                if (info != 0)
                    return linalg.Error.FactorizationFailed;

                return .init(n, ud.data, ipiv.ptr, types.isHermitianMatrix(A));
            }
        },
        else => unreachable,
    }
}

/// Convert the output a of `sytrf` or `hetrf` into the u and d matrices of the
/// UDUT factorization. On output, `a` is `u` and `d` is overwritten. `p` is
/// used as workspace for the permutation.
fn uconvert(
    comptime order: Order,
    n: u32,
    a: anytype,
    u: anytype,
    d: anytype,
    p: [*]u32,
    ipiv: [*]i32,
    hermitian: bool,
    ctx: anytype,
) !void {
    var ii: i32 = types.scast(i32, n - 1);
    while (ii >= 0) : (ii -= 1) {
        var i: u32 = types.scast(u32, ii);
        if (ipiv[i] > 0) {
            // 1x1 pivot

            // Update p
            var tmp: u32 = p[i];
            p[i] = p[types.scast(u32, ipiv[i] - 1)];
            p[types.scast(u32, ipiv[i] - 1)] = tmp;

            // Extract d
            if (comptime @TypeOf(d) != @TypeOf(null)) {
                d[i + (n - 1)] = try ops.copy(a[utils.cindex(order, i, i, n)], ctx); // diagonal

                if (i < n - 1) { // zero out under d11
                    d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                    d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
                }
            }

            // Extract u
            if (comptime @TypeOf(u) != @TypeOf(null)) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    u[utils.cindex(order, j, i, n)] = try ops.copy(a[utils.cindex(order, j, i, n)], ctx);
                }
            }

            // Permute rows of u
            while (p[i] != i) {
                const jrow: u32 = p[i];

                if (comptime @TypeOf(u) != @TypeOf(null))
                    swap_rows_suffix(order, n, u, i, jrow, i + 1);

                tmp = p[i];
                p[i] = p[jrow];
                p[jrow] = tmp;
            }
        } else {
            // 2x2 pivot

            // Update p
            var tmp: u32 = p[i - 1];
            p[i - 1] = p[types.scast(u32, -ipiv[i] - 1)];
            p[types.scast(u32, -ipiv[i] - 1)] = tmp;

            // Extract d
            if (comptime @TypeOf(d) != @TypeOf(null)) {
                d[(i - 1) + (n - 1)] = try ops.copy(a[utils.cindex(order, i - 1, i - 1, n)], ctx); // diagonal (a11)
                d[i - 1] = if (hermitian)
                    try ops.conj(a[utils.cindex(order, i - 1, i, n)], ctx)
                else
                    try ops.copy(a[utils.cindex(order, i - 1, i, n)], ctx); // subdiagonal (a21)
                d[(i - 1) + (n - 1) + n] = try ops.copy(a[utils.cindex(order, i - 1, i, n)], ctx); // superdiagonal (a12)
                d[i + (n - 1)] = try ops.copy(a[utils.cindex(order, i, i, n)], ctx); // diagonal (a22)

                if (i < n - 1) { // zero out under d22
                    d[i] = try constants.zero(Child(@TypeOf(a)), ctx); // subdiagonal
                    d[i + (n - 1) + n] = try constants.zero(Child(@TypeOf(a)), ctx); // superdiagonal
                }
            }

            // Extract u
            if (comptime @TypeOf(u) != @TypeOf(null)) {
                var j: u32 = 0;
                while (j < i - 1) : (j += 1) {
                    u[utils.cindex(order, j, i - 1, n)] = try ops.copy(a[utils.cindex(order, j, i - 1, n)], ctx);
                    u[utils.cindex(order, j, i, n)] = try ops.copy(a[utils.cindex(order, j, i, n)], ctx);
                }

                u[utils.cindex(order, i - 1, i, n)] = try constants.zero(Child(@TypeOf(a)), ctx);
            }

            ii -= 1;
            i -= 1;

            // Permute rows of u
            while (p[i] != i) {
                const jrow: u32 = p[i];

                if (comptime @TypeOf(u) != @TypeOf(null))
                    swap_rows_suffix(order, n, u, i, jrow, i + 2);

                tmp = p[i];
                p[i] = p[jrow];
                p[jrow] = tmp;
            }
        }
    }
}

// swap rows r1 <-> r2 across columns [from_col .. n-1]
fn swap_rows_suffix(comptime order: Order, n: u32, a: anytype, r1: u32, r2: u32, from_col: u32) void {
    var j: u32 = from_col;
    while (j < n) : (j += 1) {
        const tmp = a[utils.cindex(order, r1, j, n)];
        a[utils.cindex(order, r1, j, n)] = a[utils.cindex(order, r2, j, n)];
        a[utils.cindex(order, r2, j, n)] = tmp;
    }
}

fn upermute(n: u32, p: [*]u32, ipiv: [*]i32) void {
    var i: i32 = types.scast(i32, n - 1);
    while (i >= 0) : (i -= 1) {
        if (ipiv[types.scast(u32, i)] > 0) {
            // 1x1 pivot
            const tmp: u32 = p[types.scast(u32, i)];
            p[types.scast(u32, i)] = p[types.scast(u32, ipiv[types.scast(u32, i)] - 1)];
            p[types.scast(u32, ipiv[types.scast(u32, i)] - 1)] = tmp;
        } else {
            // 2x2 pivot
            const tmp: u32 = p[types.scast(u32, i - 1)];
            p[types.scast(u32, i - 1)] = p[types.scast(u32, -ipiv[types.scast(u32, i)] - 1)];
            p[types.scast(u32, -ipiv[types.scast(u32, i)] - 1)] = tmp;

            i -= 1;
        }
    }
}

pub fn bunchkaufman(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !(if (types.uploOf(@TypeOf(a)) == .lower) LDLT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) else UDUT(Numeric(@TypeOf(a)), orderOf(@TypeOf(a)))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .symmetric and types.matrixType(A) != .hermitian))
        @compileError("bunchkaufman: argument must be a symmetric or hermitian matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.bunchkaufman not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.uploOf(A)) {
        .lower => return ldlt(allocator, a, ctx),
        .upper => return udut(allocator, a, ctx),
    }
}
