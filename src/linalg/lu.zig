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

const utils = @import("utils.zig");

pub fn LU(M: type) type {
    // M is a matrix type. It conditions the parameters and information stored about the decomposition.
    return struct {
        l: matrix.Triangular(T, .lower, .unit, order),
        u: matrix.Triangular(T, .upper, .non_unit, order),

        pub fn init(lu: anytype, m: u32, n: u32) LU(Numeric(Child(@TypeOf(lu))), order) {
            return .{
                .l = .{
                    .data = lu,
                    .rows = m,
                    .cols = int.min(m, n),
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = true },
                },
                .u = .{
                    .data = lu,
                    .rows = int.min(m, n),
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = false },
                },
            };
        }

        pub fn deinit(self: *LU(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self.l.data[0 .. self.l.rows * self.u.cols]);

            self.* = undefined;
        }

        pub fn reconstruct(self: *const LU(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = L * U
            @compileError("LU.reconstruct not implemented yet");
        }
    };
}

pub fn PLU(T: type, order: Order) type {
    return struct {
        _rows: u32,
        _cols: u32,
        _ipiv: [*]i32, // 1-based indexing
        _lu: [*]T,

        pub fn init(ipiv: [*]i32, a: anytype, rows: u32, cols: u32) PLU(Numeric(Child(@TypeOf(a))), order) {
            return .{
                ._rows = rows,
                ._cols = cols,
                ._ipiv = ipiv,
                ._lu = a,
            };
        }

        pub fn deinit(self: *PLU(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._ipiv[0..int.min(self._rows, self._cols)]);
            allocator.free(self._lu[0 .. self._rows * self._cols]);

            self.* = undefined;
        }

        pub fn p(self: *const PLU(T, order), allocator: std.mem.Allocator) !matrix.Permutation(T) {
            var _p: matrix.Permutation(T) = try .init(allocator, self._rows);
            errdefer _p.deinit(allocator);

            var i: i32 = 0;
            while (i < self._rows) : (i += 1) {
                _p.data[types.scast(u32, i)] = types.scast(u32, i);
            }

            i = types.scast(i32, int.min(self._rows, self._cols) - 1);
            while (i >= 0) : (i -= 1) {
                const tmp: u32 = _p.data[types.scast(u32, i)];
                _p.data[types.scast(u32, i)] = _p.data[types.scast(u32, self._ipiv[types.scast(u32, i)] - 1)];
                _p.data[types.scast(u32, self._ipiv[types.scast(u32, i)] - 1)] = tmp;
            }

            return _p;
        }

        pub fn l(self: *const PLU(T, order)) matrix.Triangular(T, .lower, .unit, order) {
            return .{
                .data = self._lu,
                .rows = self._rows,
                .cols = int.min(self._rows, self._cols),
                .ld = if (order == .col_major) self._rows else self._cols,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn u(self: *const PLU(T, order)) matrix.Triangular(T, .upper, .non_unit, order) {
            return .{
                .data = self._lu,
                .rows = int.min(self._rows, self._cols),
                .cols = self._cols,
                .ld = if (order == .col_major) self._rows else self._cols,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn reconstruct(self: *const PLU(T, order), allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
            // A = P * L * U
            @compileError("PLU.reconstruct not implemented yet");
        }
    };
}

pub fn PLUQ(T: type, order: Order) type {
    return struct {
        p: matrix.Permutation(T),
        l: matrix.Triangular(T, .lower, .unit, order),
        u: matrix.Triangular(T, .upper, .non_unit, order),
        q: matrix.Permutation(T),

        pub fn init(p: [*]u32, lu: anytype, q: [*]u32, m: u32, n: u32) PLUQ(Numeric(Child(@TypeOf(lu))), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = m,
                    .direction = .forward,
                    .flags = .{ .owns_data = true },
                },
                .l = .{
                    .data = lu,
                    .rows = m,
                    .cols = int.min(m, n),
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = true },
                },
                .u = .{
                    .data = lu,
                    .rows = int.min(m, n),
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = false },
                },
                .q = .{
                    .data = q,
                    .size = n,
                    .direction = .forward,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *PLUQ(T, order), allocator: std.mem.Allocator) void {
            self.p.deinit(allocator);
            allocator.free(self.l.data[0 .. self.l.rows * self.u.cols]);
            self.q.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn plu(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !PLU(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or (types.matrixType(A) != .general and types.matrixType(A) != .banded and types.matrixType(A) != .tridiagonal))
        @compileError("plu: argument must be a general, banded or tridiagonal matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.plu not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .general => {
            const m: u32 = a.rows;
            const n: u32 = a.cols;

            var lu: matrix.General(Numeric(@TypeOf(a)), orderOf(A)) = try a.copy(allocator, ctx);
            errdefer lu.deinit(allocator);

            const ipiv: []i32 = try allocator.alloc(i32, int.min(m, n));
            errdefer allocator.free(ipiv);

            const info: i32 = try linalg.lapack.getrf(
                types.orderOf(A),
                types.scast(i32, m),
                types.scast(i32, n),
                lu.data,
                types.scast(i32, lu.ld),
                ipiv.ptr,
                ctx,
            );

            if (info != 0)
                return error.SingularMatrix;

            return .init(ipiv.ptr, lu.data, m, n);
        },
        .banded => return linalg.Error.NotImplemented, // lapack.gbtrf
        .tridiagonal => return linalg.Error.NotImplemented, // lapack.gttrf
        else => unreachable,
    }
}

pub fn pluq(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !PLUQ(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or types.matrixType(A) != .general)
        @compileError("pluq: argument must be a general matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.pluq not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    var lu: matrix.General(Numeric(@TypeOf(a)), orderOf(A)) = try a.copy(allocator, ctx);
    errdefer lu.deinit(allocator);

    const ipiv: []i32 = try allocator.alloc(i32, a.rows);
    errdefer allocator.free(ipiv);

    const q: []u32 = try allocator.alloc(u32, a.cols);
    errdefer allocator.free(q);

    const nb: u32 = 64; // Block size, tune for performance

    try k_pluq(
        a.rows,
        a.cols,
        nb,
        lu.data,
        lu.ld,
        ipiv.ptr,
        q.ptr,
        ctx,
    );

    var p: []u32 = try allocator.alloc(u32, a.rows);
    errdefer allocator.free(p);

    // Invert p
    var i: u32 = 0;
    while (i < a.rows) : (i += 1) {
        p.data[ipiv.data[i]] = i;
    }

    return .init(p.ptr, lu.data, q.ptr, a.rows, a.cols);
}

fn k_pluq(
    order: Order,
    m: u32,
    n: u32,
    nb: u32,
    a: anytype,
    lda: u32,
    p: [*]u32,
    q: [*]u32,
    ctx: anytype,
) !void {
    const min_dim: u32 = int.min(m, n);

    var i: u32 = 0;
    while (i < min_dim) : (i += 1) {
        p[i] = i;
        q[i] = i;
    }

    if (m < n) {
        while (i < n) : (i += 1) {
            q[i] = i;
        }
    } else if (n < m) {
        while (i < m) : (i += 1) {
            p[i] = i;
        }
    }

    if (min_dim <= nb) {
        return k_pluq2(order, m, n, a, lda, p, q, ctx);
    } else { // Remove when blocked algorithm works
        return k_pluq2(order, m, n, a, lda, p, q, ctx);
    }

    // Unreachable: blocked algorithm does not work; must be reworked

    var k: u32 = 0;
    while (k < min_dim) : (k += nb) {
        const kb: u32 = int.min(min_dim - k, nb);

        var pb: [nb]u32 = blk: {
            comptime var arr: [nb]u32 = .{0} ** nb;
            comptime var j: u32 = 0;
            inline while (j < nb) : (j += 1) {
                arr[j] = j;
            }
            break :blk arr;
        };
        var qb: [nb]u32 = blk: {
            comptime var arr: [nb]u32 = .{0} ** nb;
            comptime var j: u32 = 0;
            inline while (j < nb) : (j += 1) {
                arr[j] = j;
            }
            break :blk arr;
        };

        // Factor A[k:m, k:k + kb]
        try k_pluq2(order, m - k, kb, a + k + k * lda, lda, &pb, &qb, ctx);

        var pt: [nb]u32 = undefined;
        var qt: [nb]u32 = undefined;

        i = 0;
        while (i < kb) : (i += 1) {
            pt[i] = p[i + k];
            qt[i] = q[i + k];
        }

        i = 0;
        while (i < kb) : (i += 1) {
            p[i + k] = pt[pb[i]];
            q[i + k] = qt[qb[i]];
        }

        // Apply row permutations to A[k:k + kb, 0:k] and A[k:k + kb, k + kb:n]
        i = 0;
        while (i < kb) : (i += 1) {
            if (pb[i] == i) continue;

            // Swap A[k + i, 0:k] with A[k + pb[i], 0:k]
            if (k > 0) {
                try linalg.blas.swap(
                    types.scast(i32, k),
                    a + i + k,
                    utils.row_ld(order, lda),
                    a + pb[i] + k,
                    utils.row_ld(order, lda),
                    ctx,
                );
            }

            if (n > k + kb) {
                // Swap A[k + i, k + kb:n] with A[k + pb[i], k + kb:n]
                try linalg.blas.swap(
                    types.scast(i32, n - (k + kb)),
                    a + i + k + (k + kb) * lda,
                    utils.row_ld(order, lda),
                    a + pb[i] + k + (k + kb) * lda,
                    utils.row_ld(order, lda),
                    ctx,
                );
            }
        }

        // Apply column permutations to A[0:k, k:k + kb] and A[k + kb:m, k:k + kb]
        var j: u32 = 0;
        while (j < kb) : (j += 1) {
            if (qb[j] == j) continue;

            // Swap A[0:k, k + j] with A[0:k, k + qb[j]]
            try linalg.blas.swap(
                types.scast(i32, k),
                a + (j + k) * lda,
                utils.col_ld(order, lda),
                a + (qb[j] + k) * lda,
                utils.col_ld(order, lda),
                ctx,
            );

            // Swap A[k + kb:m, k + j] with A[k + kb:m, k + qb[j]]
            if (m > k + kb) {
                try linalg.blas.swap(
                    types.scast(i32, m - (k + kb)),
                    a + (k + kb) + (j + k) * lda,
                    utils.col_ld(order, lda),
                    a + (k + kb) + (qb[j] + k) * lda,
                    utils.col_ld(order, lda),
                    ctx,
                );
            }
        }

        // Update trailing matrix A[k + kb:m, k + kb:n]
        if (k + kb < min_dim) {
            // Solve L11 * U12 = A12, with U12 = A[k:k + kb, k + kb:n]
            if (k + kb < n) {
                try linalg.blas.trsm(
                    .col_major,
                    .left,
                    .lower,
                    .no_trans,
                    .unit,
                    types.scast(i32, kb),
                    types.scast(i32, n - (k + kb)),
                    1,
                    a + k + k * lda,
                    types.scast(i32, lda),
                    a + k + (k + kb) * lda,
                    types.scast(i32, lda),
                    ctx,
                );
            }

            // Solve L21 * U11 = A21, with L21 = A[k + kb:m, k:k + kb]
            if (k + kb < m) {
                try linalg.blas.trsm(
                    .col_major,
                    .right,
                    .upper,
                    .no_trans,
                    .non_unit,
                    types.scast(i32, m - (k + kb)),
                    types.scast(i32, kb),
                    1,
                    a + k + k * lda,
                    types.scast(i32, lda),
                    a + (k + kb) + k * lda,
                    types.scast(i32, lda),
                    ctx,
                );
            }

            if (k + kb < m and k + kb < n) {
                // Update Schur complement A22 = A22 - L21 * U12
                try linalg.blas.gemm(
                    .col_major,
                    .no_trans,
                    .no_trans,
                    types.scast(i32, m - (k + kb)),
                    types.scast(i32, n - (k + kb)),
                    types.scast(i32, kb),
                    -1,
                    a + (k + kb) + k * lda,
                    types.scast(i32, lda),
                    a + k + (k + kb) * lda,
                    types.scast(i32, lda),
                    1,
                    a + (k + kb) + (k + kb) * lda,
                    types.scast(i32, lda),
                    ctx,
                );
            }
        }
    }
}

fn k_pluq2(
    order: Order,
    m: u32,
    n: u32,
    a: anytype,
    lda: u32,
    p: [*]u32,
    q: [*]u32,
    ctx: anytype,
) !void {
    const min_dim: u32 = if (m < n) m else n;

    var k: u32 = 0;
    while (k < min_dim) : (k += 1) {
        // Find pivot in column k and swap rows
        var i_max: u32 = k;
        var j_max: u32 = k;
        var v_max: Numeric(Child(@TypeOf(a))) = ops.abs1(a[k + k * lda], ctx) catch unreachable;
        var j: u32 = k;
        while (j < n) : (j += 1) {
            const i_max_j: u32 = linalg.blas.iamax(
                types.scast(i32, m - k),
                a + k + j * lda,
                utils.col_ld(order, lda),
                ctx,
            ) catch unreachable;

            const v: Numeric(Child(@TypeOf(a))) = ops.abs1(a[k + i_max_j + j * lda], ctx) catch unreachable;
            if (ops.gt(v, v_max, ctx) catch unreachable) {
                v_max = v;
                i_max = k + i_max_j;
                j_max = j;
            }
        }

        if (ops.eq(v_max, 0, ctx) catch unreachable)
            return linalg.Error.SingularMatrix;

        // Swap rows
        if (i_max != k) {
            try linalg.blas.swap(
                types.scast(i32, n),
                a + k,
                utils.row_ld(order, lda),
                a + i_max,
                utils.row_ld(order, lda),
                ctx,
            );

            const tmp: u32 = p[k];
            p[k] = p[i_max];
            p[i_max] = tmp;
        }

        // Swap columns
        if (j_max != k) {
            try linalg.blas.swap(
                types.scast(i32, m),
                a + k * lda,
                utils.col_ld(order, lda),
                a + j_max * lda,
                utils.col_ld(order, lda),
                ctx,
            );

            const tmp: u32 = q[k];
            q[k] = q[j_max];
            q[j_max] = tmp;
        }

        // Update A[k + 1:m, k] = A[k + 1:m, k] / A[k, k]
        if (k + 1 < m) {
            const oo_akk: Numeric(Child(@TypeOf(a))) = try ops.div(1, a[k + k * lda], ctx);
            try linalg.blas.scal(
                types.scast(i32, m - (k + 1)),
                oo_akk,
                a + (k + 1) + k * lda,
                utils.col_ld(order, lda),
                ctx,
            );

            if (k + 1 < n) {
                // Update A[k + 1:m, k + 1:n] -= A[k + 1:m, k] * A[k, k + 1:n]
                try linalg.blas.ger(
                    .col_major,
                    types.scast(i32, m - (k + 1)),
                    types.scast(i32, n - (k + 1)),
                    -1,
                    a + (k + 1) + k * lda,
                    utils.col_ld(order, lda),
                    a + k + (k + 1) * lda,
                    utils.row_ld(order, lda),
                    a + (k + 1) + (k + 1) * lda,
                    types.scast(i32, lda),
                    ctx,
                );
            }
        }
    }
}
