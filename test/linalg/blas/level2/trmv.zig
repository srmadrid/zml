const std = @import("std");

const zsl = @import("zsl");
const trmv = zsl.linalg.blas.trmv;
const trmvParallel = zsl.linalg.blas.trmvParallel;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 0), @as(usize, 1), @as(isize, 1) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 0), @as(usize, 1), @as(isize, -1) },

    // Single element
    .{ .row_major, .upper, .conj_no_trans, .unit, @as(usize, 1), @as(usize, 1), @as(isize, 1) },
    .{ .col_major, .lower, .no_trans, .non_unit, @as(usize, 1), @as(usize, 1), @as(isize, -1) },

    // Square aligned
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 32), @as(usize, 32), @as(isize, 1) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 64), @as(usize, 64), @as(isize, 1) },

    // Unaligned sizes
    .{ .col_major, .lower, .conj_trans, .non_unit, @as(usize, 37), @as(usize, 37), @as(isize, 1) },
    .{ .row_major, .upper, .no_trans, .unit, @as(usize, 69), @as(usize, 69), @as(isize, 1) },

    // Leading dimension variations (lda > n)
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 32), @as(usize, 40), @as(isize, 1) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 37), @as(usize, 45), @as(isize, 1) },

    // Strided vectors
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 32), @as(usize, 32), @as(isize, 2) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 32), @as(usize, 32), @as(isize, -3) },
    .{ .col_major, .lower, .conj_trans, .non_unit, @as(usize, 37), @as(usize, 40), @as(isize, -1) },
    .{ .row_major, .upper, .conj_no_trans, .unit, @as(usize, 37), @as(usize, 40), @as(isize, 2) },
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 64), @as(usize, 64), @as(isize, -1) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 65), @as(usize, 65), @as(isize, 7) },

    // Large matrices for parallel chunking
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 2000), @as(usize, 2000), @as(isize, 1) },
    .{ .row_major, .lower, .trans, .unit, @as(usize, 1500), @as(usize, 1500), @as(isize, 2) },
    .{ .row_major, .upper, .conj_no_trans, .non_unit, @as(usize, 1000), @as(usize, 1000), @as(isize, 2) },
    .{ .col_major, .lower, .trans, .unit, @as(usize, 1500), @as(usize, 1505), @as(isize, 1) },
    .{ .col_major, .upper, .no_trans, .non_unit, @as(usize, 1000), @as(usize, 1000), @as(isize, 1) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .lower, .no_trans, .unit, @as(usize, 1003), @as(usize, 1510), @as(isize, 1) },
};

test trmv {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeTrmvTest(
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
                &pool,
            );
        }
    }
}

fn executeTrmvTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    uplo: zsl.matrix.Uplo,
    transa: zsl.linalg.blas.Transpose,
    diag: zsl.matrix.Diag,
    n: usize,
    lda: usize,
    incx: isize,
    pool: *zsl.thread.Pool,
) !void {
    if (N == f64 and (transa == .conj_no_trans or transa == .conj_trans))
        return;

    const abs_incx = @abs(incx);
    const len_mem_x = if (n == 0) 0 else 1 + (n - 1) * abs_incx;
    const len_mem_a = lda * n;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const x_expected = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x_expected);

    const x_actual = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x_actual);
    const x_actual_parallel = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x_actual_parallel);

    const work = try allocator.alloc(N, len_mem_x);
    defer allocator.free(work);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (x_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        x_actual[i] = elem.*;
        x_actual_parallel[i] = elem.*;
    }

    if (N == f64) {
        zsl.linalg.cblas.dtrmv(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            diag.toInt(c_int),
            zsl.numeric.cast(isize, n),
            a.ptr,
            zsl.numeric.cast(isize, lda),
            x_expected.ptr,
            zsl.numeric.cast(isize, incx),
        );
    } else {
        zsl.linalg.cblas.ztrmv(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            diag.toInt(c_int),
            zsl.numeric.cast(isize, n),
            a.ptr,
            zsl.numeric.cast(isize, lda),
            x_expected.ptr,
            zsl.numeric.cast(isize, incx),
        );
    }

    trmv(
        layout,
        uplo,
        transa,
        diag,
        n,
        a.ptr,
        lda,
        x_actual.ptr,
        incx,
    ) catch |e| {
        std.debug.print("\n\tTRMV Execution Failure\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | transa: {s} | diag: {s} | n: {} | incx: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), @tagName(transa), @tagName(diag), n, incx });
        return e;
    };

    trmvParallel(
        layout,
        uplo,
        transa,
        diag,
        n,
        a.ptr,
        lda,
        x_actual_parallel.ptr,
        incx,
        work.ptr,
        pool,
    ) catch |e| {
        std.debug.print("\n\tTRMV Execution Failure\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | transa: {s} | diag: {s} | n: {} | incx: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), @tagName(transa), @tagName(diag), n, incx });
        return e;
    };

    if (n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    for (x_expected, x_actual, x_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            const diff_parallel = zsl.float.abs(expected - actual_parallel);
            if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
                (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
            {
                std.debug.print("\n\tTRMV Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | uplo: {s} | transa: {s} | diag: {s} | n: {} | incx: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transa), @tagName(diag), n, incx, i });
                std.debug.print("Expected (CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (zsl):       {d}\n", .{actual});
                std.debug.print("Diff:               {d}\n", .{diff});
                std.debug.print("Actual (zsl, pool): {d}\n", .{actual_parallel});
                std.debug.print("Diff (pool):        {d}\n", .{diff_parallel});

                return error.TestFailed;
            }
        } else {
            const diff_re = zsl.float.abs(expected.re - actual.re);
            const diff_im = zsl.float.abs(expected.im - actual.im);
            const diff_parallel_re = zsl.float.abs(expected.re - actual_parallel.re);
            const diff_parallel_im = zsl.float.abs(expected.im - actual_parallel.im);

            const mag_re = zsl.float.abs(expected.re);
            const mag_im = zsl.float.abs(expected.im);

            if ((diff_re > abs_tol and diff_re > mag_re * rel_tol) or
                (diff_im > abs_tol and diff_im > mag_im * rel_tol) or
                (diff_parallel_re > abs_tol and diff_parallel_re > mag_re * rel_tol) or
                (diff_parallel_im > abs_tol and diff_parallel_im > mag_im * rel_tol))
            {
                std.debug.print("\n\tTRMV Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | uplo: {s} | transa: {s} | diag: {s} | n: {} | incx: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transa), @tagName(diag), n, incx, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });
                std.debug.print("Actual (zsl, pool): {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff (pool):        {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
