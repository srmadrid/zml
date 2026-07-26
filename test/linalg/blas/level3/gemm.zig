const std = @import("std");

const zsl = @import("zsl");
const gemm = zsl.linalg.blas.gemm;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays (all zero)
    .{ .col_major, .no_trans, .no_trans, @as(usize, 0), @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .trans, .trans, @as(usize, 0), @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // k = 0
    .{ .col_major, .no_trans, .no_trans, @as(usize, 10), @as(usize, 15), @as(usize, 0), @as(usize, 10), @as(usize, 1), @as(usize, 10) },
    .{ .row_major, .trans, .no_trans, @as(usize, 12), @as(usize, 8), @as(usize, 0), @as(usize, 12), @as(usize, 8), @as(usize, 8) },

    // m = 0 or n = 0
    .{ .col_major, .no_trans, .trans, @as(usize, 0), @as(usize, 10), @as(usize, 5), @as(usize, 1), @as(usize, 10), @as(usize, 1) },
    .{ .row_major, .conj_trans, .no_trans, @as(usize, 10), @as(usize, 0), @as(usize, 5), @as(usize, 10), @as(usize, 1), @as(usize, 1) },

    // Single element
    .{ .col_major, .no_trans, .no_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .conj_no_trans, .conj_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .no_trans, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32) },
    .{ .row_major, .trans, .trans, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64) },

    // Unaligned / Rectangular shapes
    .{ .col_major, .no_trans, .no_trans, @as(usize, 13), @as(usize, 17), @as(usize, 19), @as(usize, 13), @as(usize, 19), @as(usize, 13) },
    .{ .row_major, .no_trans, .no_trans, @as(usize, 17), @as(usize, 19), @as(usize, 13), @as(usize, 13), @as(usize, 19), @as(usize, 19) },
    .{ .col_major, .trans, .no_trans, @as(usize, 13), @as(usize, 17), @as(usize, 19), @as(usize, 19), @as(usize, 19), @as(usize, 13) },
    .{ .row_major, .no_trans, .trans, @as(usize, 13), @as(usize, 17), @as(usize, 19), @as(usize, 19), @as(usize, 19), @as(usize, 17) },

    // Complex Conjugation variations
    .{ .col_major, .conj_trans, .no_trans, @as(usize, 15), @as(usize, 12), @as(usize, 10), @as(usize, 10), @as(usize, 10), @as(usize, 15) },
    .{ .col_major, .no_trans, .conj_trans, @as(usize, 12), @as(usize, 15), @as(usize, 10), @as(usize, 12), @as(usize, 15), @as(usize, 12) },
    .{ .row_major, .conj_no_trans, .conj_no_trans, @as(usize, 14), @as(usize, 16), @as(usize, 18), @as(usize, 18), @as(usize, 16), @as(usize, 16) },
    .{ .row_major, .conj_trans, .conj_trans, @as(usize, 16), @as(usize, 14), @as(usize, 18), @as(usize, 16), @as(usize, 18), @as(usize, 14) },

    // Padded leading dimensions (lda, ldb, ldc > minimum required)
    .{ .col_major, .no_trans, .no_trans, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 16), @as(usize, 18), @as(usize, 20) },
    .{ .row_major, .trans, .no_trans, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 15), @as(usize, 15), @as(usize, 16) },
    .{ .col_major, .conj_trans, .conj_trans, @as(usize, 8), @as(usize, 9), @as(usize, 10), @as(usize, 12), @as(usize, 12), @as(usize, 12) },

    // Thin / Tall matrices
    .{ .col_major, .no_trans, .no_trans, @as(usize, 500), @as(usize, 2), @as(usize, 10), @as(usize, 500), @as(usize, 10), @as(usize, 500) },
    .{ .row_major, .no_trans, .no_trans, @as(usize, 2), @as(usize, 500), @as(usize, 10), @as(usize, 10), @as(usize, 500), @as(usize, 500) },
    .{ .col_major, .trans, .trans, @as(usize, 10), @as(usize, 10), @as(usize, 500), @as(usize, 500), @as(usize, 10), @as(usize, 10) },

    // Large arrays
    .{ .col_major, .no_trans, .no_trans, @as(usize, 120), @as(usize, 100), @as(usize, 150), @as(usize, 120), @as(usize, 150), @as(usize, 120) },
    .{ .row_major, .trans, .no_trans, @as(usize, 100), @as(usize, 120), @as(usize, 150), @as(usize, 150), @as(usize, 150), @as(usize, 120) },
    .{ .col_major, .no_trans, .trans, @as(usize, 120), @as(usize, 150), @as(usize, 100), @as(usize, 120), @as(usize, 150), @as(usize, 120) },
    .{ .row_major, .trans, .trans, @as(usize, 150), @as(usize, 100), @as(usize, 120), @as(usize, 150), @as(usize, 120), @as(usize, 100) },
    .{ .row_major, .trans, .trans, @as(usize, 1500), @as(usize, 1000), @as(usize, 1200), @as(usize, 1500), @as(usize, 1200), @as(usize, 1000) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .no_trans, .no_trans, @as(usize, 67), @as(usize, 71), @as(usize, 73), @as(usize, 70), @as(usize, 75), @as(usize, 70) },
    .{ .row_major, .trans, .conj_trans, @as(usize, 71), @as(usize, 67), @as(usize, 73), @as(usize, 75), @as(usize, 75), @as(usize, 67) },
};

test gemm {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeGemmTest(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                combo[3],
                combo[4],
                combo[5],
                combo[6],
                combo[7],
                combo[8],
                &pool,
            );
        }
    }
}

fn executeGemmTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    transa: zsl.linalg.blas.Transpose,
    transb: zsl.linalg.blas.Transpose,
    m: usize,
    n: usize,
    k: usize,
    lda: usize,
    ldb: usize,
    ldc: usize,
    pool: *zsl.thread.Pool,
) !void {
    _ = pool;

    if (N == f64 and (transa == .conj_no_trans or transa == .conj_trans or
        transb == .conj_no_trans or transb == .conj_trans))
        return;

    const notransa = transa == .no_trans or transa == .conj_no_trans;
    const notransb = transb == .no_trans or transb == .conj_no_trans;

    const len_mem_a = if (layout == .col_major)
        if (notransa) lda * k else lda * m
    else if (notransa) lda * m else lda * k;

    const len_mem_b = if (layout == .col_major)
        if (notransb) ldb * n else ldb * k
    else if (notransb) ldb * k else ldb * n;

    const len_mem_c = if (layout == .col_major) ldc * n else ldc * m;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const b = try allocator.alloc(N, len_mem_b);
    defer allocator.free(b);

    const c_expected = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_expected);

    const c_actual = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_actual);
    const c_actual_parallel = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_actual_parallel);

    const alpha: N = tzsl.randomNumber(N, rand);
    const beta: N = tzsl.randomNumber(N, rand);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (b) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (c_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        c_actual[i] = elem.*;
        c_actual_parallel[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dgemm(
            layout.toInt(c_int),
            transa.toInt(c_int),
            transb.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b.ptr,
            zsl.numeric.cast(isize, ldb),
            beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        )
    else
        zsl.linalg.cblas.zgemm(
            layout.toInt(c_int),
            transa.toInt(c_int),
            transb.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b.ptr,
            zsl.numeric.cast(isize, ldb),
            &beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );

    gemm(
        layout,
        transa,
        transb,
        m,
        n,
        k,
        alpha,
        a.ptr,
        lda,
        b.ptr,
        ldb,
        beta,
        c_actual.ptr,
        ldc,
    ) catch |e| {
        std.debug.print("\n\tGEMM Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | transa: {s} | transb: {s} | m: {} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(transa), @tagName(transb), m, n, k, lda, ldb, ldc });
        return e;
    };

    // gemmParallel(
    //     layout,
    //     transa,
    //     transb,
    //     m,
    //     n,
    //     k,
    //     alpha,
    //     a.ptr,
    //     lda,
    //     b.ptr,
    //     ldb,
    //     beta,
    //     c_actual_parallel.ptr,
    //     ldc,
    //     pool,
    // ) catch |e| {
    //     std.debug.print("\n\tGEMM Test Failed\n", .{});
    //     std.debug.print("Type: {s} | layout: {s} | transa: {s} | transb: {s} | m: {} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(transa), @tagName(transb), m, n, k, lda, ldb, ldc });
    //     return e;
    // };

    if (m == 0 or n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, zsl.int.max(1, k));
    const abs_tol = 1e-10;
    for (c_expected, c_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tGEMM Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | transa: {s} | transb: {s} | m: {} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(transa), @tagName(transb), m, n, k, lda, ldb, ldc, i });
                std.debug.print("Expected (CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (zsl):       {d}\n", .{actual});
                std.debug.print("Diff:               {d}\n", .{diff});

                return error.TestFailed;
            }
        } else {
            const diff_re = zsl.float.abs(expected.re - actual.re);
            const diff_im = zsl.float.abs(expected.im - actual.im);

            const mag_re = zsl.float.abs(expected.re);
            const mag_im = zsl.float.abs(expected.im);

            if ((diff_re > abs_tol and diff_re > mag_re * rel_tol) or
                (diff_im > abs_tol and diff_im > mag_im * rel_tol))
            {
                std.debug.print("\n\tGEMM Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | transa: {s} | transb: {s} | m: {} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(transa), @tagName(transb), m, n, k, lda, ldb, ldc, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
