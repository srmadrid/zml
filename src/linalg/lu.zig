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
const General = matrix.General;
const Triangular = matrix.Triangular;
const Permutation = matrix.Permutation;

const linalg = @import("../linalg.zig");

pub fn LU(T: type, order: Order) type {
    return struct {
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),

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
    };
}

pub fn PLU(T: type, order: Order) type {
    return struct {
        p: Permutation(T),
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),

        pub fn init(p: [*]u32, lu: anytype, m: u32, n: u32) PLU(Numeric(Child(@TypeOf(lu))), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = m,
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
            };
        }

        pub fn deinit(self: *PLU(T, order), allocator: std.mem.Allocator) void {
            self.p.deinit(allocator);
            allocator.free(self.l.data[0 .. self.l.rows * self.u.cols]);

            self.* = undefined;
        }
    };
}

pub fn PLUQ(T: type, order: Order) type {
    return struct {
        p: Permutation(T),
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),
        q: Permutation(T),

        pub fn init(p: [*]u32, lu: anytype, q: [*]u32, m: u32, n: u32) PLUQ(Numeric(Child(@TypeOf(lu))), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = m,
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

    comptime if (!types.isMatrix(A))
        @compileError("plu: argument must be a matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.matmul not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    switch (comptime types.matrixType(A)) {
        .general => {
            const m: u32 = a.rows;
            const n: u32 = a.cols;

            var lu: General(Numeric(@TypeOf(a)), orderOf(A)) = try .init(allocator, m, n);
            errdefer lu.deinit(allocator);

            if (comptime types.orderOf(A) == .col_major) {
                var j: u32 = 0;
                while (j < n) : (j += 1) {
                    try linalg.blas.copy(
                        types.scast(i32, m),
                        a.data + j * a.ld,
                        1,
                        lu.data + j * lu.ld,
                        1,
                        ctx,
                    );
                }
            } else {
                var i: u32 = 0;
                while (i < m) : (i += 1) {
                    try linalg.blas.copy(
                        types.scast(i32, n),
                        a.data + i * a.ld,
                        1,
                        lu.data + i * lu.ld,
                        1,
                        ctx,
                    );
                }
            }

            var ipiv: vector.Vector(i32) = try vector.Vector(i32).init(allocator, m);
            defer ipiv.deinit(allocator);

            const info: i32 = try linalg.lapack.getrf(
                types.orderOf(A),
                types.scast(i32, m),
                types.scast(i32, n),
                lu.data,
                types.scast(i32, lu.ld),
                ipiv.data,
                ctx,
            );

            if (info != 0)
                return error.SingularMatrix;

            // Convert from 1 to 0-based indexing and invert permutation to be
            // able to reconstruct A as P * L * U instead of P^T * L * U
            var p: vector.Vector(u32) = try vector.Vector(u32).init(allocator, m);
            errdefer p.deinit(allocator);

            var i: u32 = 0;
            while (i < m) : (i += 1) {
                p.data[i] = i;
            }

            i = m - 1;
            while (i > 0) : (i -= 1) {
                const tmp: u32 = p.data[i];
                p.data[i] = p.data[types.scast(u32, ipiv.data[i] - 1)];
                p.data[types.scast(u32, ipiv.data[i] - 1)] = tmp;
            }

            // Extra i = 0 step
            const tmp: u32 = p.data[0];
            p.data[0] = p.data[types.scast(u32, ipiv.data[0] - 1)];
            p.data[types.scast(u32, ipiv.data[0] - 1)] = tmp;

            return .init(p.data, lu.data, m, n);
        },
        .banded => return linalg.Error.NotImplemented, // lapack.gbtrf
        .tridiagonal => return linalg.Error.NotImplemented, // lapack.gttrf
        else => @compileError("plu: argument must be a general, banded or tridiagonal matrix, got " ++ @typeName(A)),
    }
}

pub fn pluq(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !PLUQ(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or types.matrixType(A) != .general)
        @compileError("pluq: argument must be a general matrix, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.matmul not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    var lu: General(Numeric(@TypeOf(a)), orderOf(A)) = try .init(allocator, a.rows, a.cols);
    errdefer lu.deinit(allocator);

    if (comptime types.orderOf(A) == .col_major) {
        var j: u32 = 0;
        while (j < a.cols) : (j += 1) {
            try linalg.blas.copy(
                types.scast(i32, a.rows),
                a.data + j * a.ld,
                1,
                lu.data + j * lu.ld,
                1,
                ctx,
            );
        }
    } else {
        var i: u32 = 0;
        while (i < a.rows) : (i += 1) {
            try linalg.blas.copy(
                types.scast(i32, a.cols),
                a.data + i * a.ld,
                1,
                lu.data + i * lu.ld,
                1,
                ctx,
            );
        }
    }

    var ipiv: vector.Vector(u32) = try vector.Vector(u32).init(allocator, a.rows);
    defer ipiv.deinit(allocator);
    var q: vector.Vector(u32) = try vector.Vector(u32).init(allocator, a.cols);
    errdefer q.deinit(allocator);

    const nb: u32 = 64; // Block size, tune for performance

    if (comptime types.orderOf(A) == .col_major) {
        try k_pluq_c(
            a.rows,
            a.cols,
            nb,
            lu.data,
            lu.ld,
            ipiv.data,
            q.data,
            ctx,
        );
    } else {
        return linalg.Error.NotImplemented;
    }

    var p: vector.Vector(u32) = try vector.Vector(u32).init(allocator, a.rows);
    errdefer p.deinit(allocator);

    // Invert p
    var i: u32 = 0;
    while (i < a.rows) : (i += 1) {
        p.data[ipiv.data[i]] = i;
    }

    return .init(p.data, lu.data, q.data, a.rows, a.cols);
}

fn k_pluq_c(
    m: u32,
    n: u32,
    comptime nb: u32,
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
        return k_pluq2_c(m, n, a, lda, p, q, ctx);
    } else {
        return k_pluq2_c(m, n, a, lda, p, q, ctx);
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
        try k_pluq2_c(m - k, kb, a + k + k * lda, lda, &pb, &qb, ctx);

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
                    types.scast(i32, lda),
                    a + pb[i] + k,
                    types.scast(i32, lda),
                    ctx,
                );
            }

            if (n > k + kb) {
                // Swap A[k + i, k + kb:n] with A[k + pb[i], k + kb:n]
                try linalg.blas.swap(
                    types.scast(i32, n - (k + kb)),
                    a + i + k + (k + kb) * lda,
                    types.scast(i32, lda),
                    a + pb[i] + k + (k + kb) * lda,
                    types.scast(i32, lda),
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
                1,
                a + (qb[j] + k) * lda,
                1,
                ctx,
            );

            // Swap A[k + kb:m, k + j] with A[k + kb:m, k + qb[j]]
            if (m > k + kb) {
                try linalg.blas.swap(
                    types.scast(i32, m - (k + kb)),
                    a + (k + kb) + (j + k) * lda,
                    1,
                    a + (k + kb) + (qb[j] + k) * lda,
                    1,
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

fn k_pluq2_c(
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
        var v_max: Numeric(Child(@TypeOf(a))) = ops.abs(a[k + k * lda], ctx) catch unreachable;
        var j: u32 = k;
        while (j < n) : (j += 1) {
            const i_max_j: u32 = linalg.blas.iamax(
                types.scast(i32, m - k),
                a + k + j * lda,
                1,
                ctx,
            ) catch unreachable;

            const v: Numeric(Child(@TypeOf(a))) = ops.abs(a[k + i_max_j + j * lda], ctx) catch unreachable;
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
                types.scast(i32, lda),
                a + i_max,
                types.scast(i32, lda),
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
                1,
                a + j_max * lda,
                1,
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
                1,
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
                    1,
                    a + k + (k + 1) * lda,
                    types.scast(i32, lda),
                    a + (k + 1) + (k + 1) * lda,
                    types.scast(i32, lda),
                    ctx,
                );
            }
        }
    }
}
