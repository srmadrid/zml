const std = @import("std");

const zsl = @import("zsl");
const trmm = zsl.linalg.blas.trmm;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays (all zero)
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .lower, .trans, .unit, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1) },

    // m = 0 or n = 0
    .{ .col_major, .left, .lower, .trans, .non_unit, @as(usize, 0), @as(usize, 10), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .upper, .no_trans, .unit, @as(usize, 10), @as(usize, 0), @as(usize, 10), @as(usize, 10) },

    // Single element
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .right, .lower, .conj_trans, .unit, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32) },
    .{ .row_major, .right, .lower, .trans, .unit, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64) },

    // Unaligned / Rectangular shapes (Left side: A is m x m)
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 13), @as(usize, 17), @as(usize, 13), @as(usize, 13) },
    .{ .row_major, .left, .lower, .no_trans, .unit, @as(usize, 17), @as(usize, 19), @as(usize, 17), @as(usize, 19) },
    .{ .col_major, .left, .lower, .no_trans, .non_unit, @as(usize, 13), @as(usize, 17), @as(usize, 19), @as(usize, 13) },

    // Unaligned / Rectangular shapes (Right side: A is n x n)
    .{ .col_major, .right, .upper, .no_trans, .non_unit, @as(usize, 13), @as(usize, 17), @as(usize, 17), @as(usize, 13) },
    .{ .row_major, .right, .lower, .trans, .unit, @as(usize, 17), @as(usize, 13), @as(usize, 13), @as(usize, 13) },

    // Complex Conjugation variations
    .{ .col_major, .left, .upper, .conj_trans, .non_unit, @as(usize, 15), @as(usize, 12), @as(usize, 15), @as(usize, 15) },
    .{ .col_major, .right, .lower, .conj_trans, .unit, @as(usize, 12), @as(usize, 15), @as(usize, 15), @as(usize, 12) },
    .{ .row_major, .left, .upper, .conj_no_trans, .non_unit, @as(usize, 14), @as(usize, 16), @as(usize, 14), @as(usize, 16) },
    .{ .row_major, .right, .lower, .conj_trans, .unit, @as(usize, 16), @as(usize, 14), @as(usize, 14), @as(usize, 14) },

    // Padded leading dimensions (lda, ldb > minimum required)
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 16) },
    .{ .row_major, .right, .lower, .trans, .unit, @as(usize, 10), @as(usize, 12), @as(usize, 15), @as(usize, 16) },
    .{ .col_major, .left, .lower, .conj_trans, .non_unit, @as(usize, 8), @as(usize, 9), @as(usize, 10), @as(usize, 12) },

    // Thin / Tall matrices
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 500), @as(usize, 2), @as(usize, 500), @as(usize, 500) },
    .{ .row_major, .right, .lower, .no_trans, .unit, @as(usize, 2), @as(usize, 500), @as(usize, 500), @as(usize, 500) },
    .{ .col_major, .left, .upper, .conj_trans, .unit, @as(usize, 10), @as(usize, 500), @as(usize, 10), @as(usize, 10) },

    // Large arrays
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 120), @as(usize, 100), @as(usize, 120), @as(usize, 120) },
    .{ .row_major, .right, .lower, .trans, .unit, @as(usize, 100), @as(usize, 120), @as(usize, 120), @as(usize, 120) },
    .{ .col_major, .left, .lower, .trans, .non_unit, @as(usize, 120), @as(usize, 150), @as(usize, 120), @as(usize, 120) },
    .{ .row_major, .right, .upper, .trans, .unit, @as(usize, 1500), @as(usize, 1000), @as(usize, 1000), @as(usize, 1000) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .left, .upper, .no_trans, .non_unit, @as(usize, 67), @as(usize, 71), @as(usize, 70), @as(usize, 70) },
    .{ .row_major, .right, .lower, .conj_trans, .unit, @as(usize, 71), @as(usize, 67), @as(usize, 67), @as(usize, 67) },
};

test trmm {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeTrmmTest(
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

fn executeTrmmTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    side: zsl.linalg.blas.Side,
    uplo: zsl.matrix.Uplo,
    transa: zsl.linalg.blas.Transpose,
    diag: zsl.matrix.Diag,
    m: usize,
    n: usize,
    lda: usize,
    ldb: usize,
    pool: *zsl.thread.Pool,
) !void {
    _ = pool;

    if (N == f64 and (transa == .conj_no_trans or transa == .conj_trans))
        return;

    const nrowa = if (side == .left) m else n;

    const len_mem_a = if (layout == .col_major) lda * nrowa else lda * nrowa;
    const len_mem_b = if (layout == .col_major) ldb * n else ldb * m;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const b_expected = try allocator.alloc(N, len_mem_b);
    defer allocator.free(b_expected);

    const b_actual = try allocator.alloc(N, len_mem_b);
    defer allocator.free(b_actual);

    const alpha: N = tzsl.randomNumber(N, rand);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (0..nrowa) |i| {
        for (0..nrowa) |j| {
            const idx = if (layout == .col_major) j * lda + i else i * lda + j;
            if (i == j and diag == .unit) {
                if (N == f64) {
                    a[idx] = 1.0;
                } else {
                    a[idx] = .{ .re = 1.0, .im = 0.0 };
                }
            } else if ((uplo == .upper and i > j) or (uplo == .lower and i < j)) {
                if (N == f64) {
                    a[idx] = 0.0;
                } else {
                    a[idx] = .{ .re = 0.0, .im = 0.0 };
                }
            }
        }
    }

    for (b_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        b_actual[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dtrmm(
            layout.toInt(c_int),
            side.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            diag.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b_expected.ptr,
            zsl.numeric.cast(isize, ldb),
        )
    else
        zsl.linalg.cblas.ztrmm(
            layout.toInt(c_int),
            side.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            diag.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b_expected.ptr,
            zsl.numeric.cast(isize, ldb),
        );

    trmm(
        layout,
        side,
        uplo,
        transa,
        diag,
        m,
        n,
        alpha,
        a.ptr,
        lda,
        b_actual.ptr,
        ldb,
    ) catch |e| {
        std.debug.print("\n\tTRMM Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | side: {s} | uplo: {s} | transa: {s} | diag: {s} | m: {} | n: {} | lda: {} | ldb: {}\n", .{ @typeName(N), @tagName(layout), @tagName(side), @tagName(uplo), @tagName(transa), @tagName(diag), m, n, lda, ldb });
        return e;
    };

    if (m == 0 or n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, zsl.int.max(1, nrowa));
    const abs_tol = 1e-10;
    for (b_expected, b_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tTRMM Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | side: {s} | uplo: {s} | transa: {s} | diag: {s} | m: {} | n: {} | lda: {} | ldb: {} | index: {}\n", .{ @tagName(layout), @tagName(side), @tagName(uplo), @tagName(transa), @tagName(diag), m, n, lda, ldb, i });
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
                std.debug.print("\n\tTRMM Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | side: {s} | uplo: {s} | transa: {s} | diag: {s} | m: {} | n: {} | lda: {} | ldb: {} | index: {}\n", .{ @tagName(layout), @tagName(side), @tagName(uplo), @tagName(transa), @tagName(diag), m, n, lda, ldb, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
