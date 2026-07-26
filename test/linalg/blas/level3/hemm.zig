const std = @import("std");

const zsl = @import("zsl");
const hemm = zsl.linalg.blas.hemm;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays (all zero)
    .{ .col_major, .left, .upper, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .lower, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // m = 0 or n = 0
    .{ .col_major, .left, .lower, @as(usize, 0), @as(usize, 10), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .upper, @as(usize, 10), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Single element
    .{ .col_major, .left, .upper, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .lower, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .left, .lower, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32) },
    .{ .row_major, .right, .upper, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64) },

    // Rectangular Left Side (A is m-by-m -> lda >= m)
    .{ .col_major, .left, .upper, @as(usize, 13), @as(usize, 19), @as(usize, 13), @as(usize, 13), @as(usize, 13) },
    .{ .row_major, .left, .lower, @as(usize, 17), @as(usize, 13), @as(usize, 17), @as(usize, 13), @as(usize, 13) },

    // Rectangular Right Side (A is n-by-n -> lda >= n)
    .{ .col_major, .right, .lower, @as(usize, 19), @as(usize, 13), @as(usize, 13), @as(usize, 19), @as(usize, 19) },
    .{ .row_major, .right, .upper, @as(usize, 13), @as(usize, 17), @as(usize, 17), @as(usize, 17), @as(usize, 17) },

    // Padded leading dimensions (lda, ldb, ldc > minimum required)
    .{ .col_major, .left, .upper, @as(usize, 10), @as(usize, 12), @as(usize, 12), @as(usize, 14), @as(usize, 15) },
    .{ .row_major, .right, .lower, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 15), @as(usize, 16) },
    .{ .col_major, .right, .upper, @as(usize, 8), @as(usize, 9), @as(usize, 12), @as(usize, 10), @as(usize, 10) },

    // Thin / Tall matrices
    .{ .col_major, .left, .lower, @as(usize, 100), @as(usize, 2), @as(usize, 100), @as(usize, 100), @as(usize, 100) },
    .{ .row_major, .right, .upper, @as(usize, 2), @as(usize, 100), @as(usize, 100), @as(usize, 100), @as(usize, 100) },

    // Large arrays
    .{ .col_major, .left, .lower, @as(usize, 120), @as(usize, 100), @as(usize, 120), @as(usize, 120), @as(usize, 120) },
    .{ .row_major, .right, .upper, @as(usize, 100), @as(usize, 120), @as(usize, 120), @as(usize, 120), @as(usize, 120) },
    .{ .row_major, .right, .upper, @as(usize, 1000), @as(usize, 1200), @as(usize, 1200), @as(usize, 1200), @as(usize, 1200) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .left, .upper, @as(usize, 67), @as(usize, 71), @as(usize, 70), @as(usize, 75), @as(usize, 70) },
    .{ .row_major, .right, .lower, @as(usize, 71), @as(usize, 67), @as(usize, 70), @as(usize, 70), @as(usize, 75) },
};

test hemm {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeHemmTest(
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
                &pool,
            );
        }
    }
}

fn executeHemmTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    side: zsl.linalg.blas.Side,
    uplo: zsl.matrix.Uplo,
    m: usize,
    n: usize,
    lda: usize,
    ldb: usize,
    ldc: usize,
    pool: *zsl.thread.Pool,
) !void {
    _ = pool;

    const nrowa: usize = if (side == .left) m else n;

    // A is a square matrix of size nrowa-by-nrowa
    const len_mem_a = lda * nrowa;

    const len_mem_b = if (layout == .col_major) ldb * n else ldb * m;
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

    if (N != f64) {
        for (0..nrowa) |i| {
            a[i * (lda + 1)].im = 0.0;
        }
    }

    for (b) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (c_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        c_actual[i] = elem.*;
        c_actual_parallel[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dsymm(
            layout.toInt(c_int),
            side.toInt(c_int),
            uplo.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
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
        zsl.linalg.cblas.zhemm(
            layout.toInt(c_int),
            side.toInt(c_int),
            uplo.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b.ptr,
            zsl.numeric.cast(isize, ldb),
            &beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );

    hemm(
        layout,
        side,
        uplo,
        m,
        n,
        alpha,
        a.ptr,
        lda,
        b.ptr,
        ldb,
        beta,
        c_actual.ptr,
        ldc,
    ) catch |e| {
        std.debug.print("\n\tHEMM Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | side: {s} | uplo: {s} | m: {} | n: {} | lda: {} | ldb: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(side), @tagName(uplo), m, n, lda, ldb, ldc });
        return e;
    };

    // hemmParallel(
    //     layout,
    //     side,
    //     uplo,
    //     m,
    //     n,
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
    //     std.debug.print("\n\tHEMM Test Failed\n", .{});
    //     std.debug.print("Type: {s} | layout: {s} | side: {s} | uplo: {s} | m: {} | n: {} | lda: {} | ldb: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(side), @tagName(uplo), m, n, lda, ldb, ldc });
    //     return e;
    // };

    if (m == 0 or n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, zsl.int.max(1, nrowa));
    const abs_tol = 1e-10;
    for (c_expected, c_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tHEMM Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | side: {s} | uplo: {s} | m: {} | n: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(side), @tagName(uplo), m, n, lda, ldb, ldc, i });
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
                std.debug.print("\n\tHEMM Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | side: {s} | uplo: {s} | m: {} | n: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(side), @tagName(uplo), m, n, lda, ldb, ldc, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
