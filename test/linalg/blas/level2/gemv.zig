const std = @import("std");

const zsl = @import("zsl");
const gemv = zsl.linalg.blas.gemv;
const gemvParallel = zsl.linalg.blas.gemvParallel;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, .no_trans, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 0), @as(usize, 10), @as(usize, 10), @as(isize, -1), @as(isize, 2) },
    .{ .col_major, .conj_trans, @as(usize, 10), @as(usize, 0), @as(usize, 10), @as(isize, 1), @as(isize, -1) },

    // Single element
    .{ .row_major, .conj_no_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(isize, 1), @as(isize, 1) },

    // Square aligned
    .{ .col_major, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(isize, 1), @as(isize, 1) },

    // Unaligned tall
    .{ .col_major, .conj_trans, @as(usize, 69), @as(usize, 37), @as(usize, 69), @as(isize, 1), @as(isize, 1) },

    // Unaligned wide
    .{ .row_major, .no_trans, @as(usize, 37), @as(usize, 69), @as(usize, 69), @as(isize, 1), @as(isize, 1) },

    // Col major requires lda >= m
    .{ .col_major, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 40), @as(isize, 1), @as(isize, 1) },
    .{ .col_major, .trans, @as(usize, 37), @as(usize, 16), @as(usize, 45), @as(isize, 1), @as(isize, 1) },

    // Row major requires lda >= n
    .{ .row_major, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 40), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .conj_no_trans, @as(usize, 16), @as(usize, 37), @as(usize, 45), @as(isize, 1), @as(isize, 1) },

    // Strided
    .{ .col_major, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(isize, 2), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 3) },
    .{ .col_major, .conj_trans, @as(usize, 37), @as(usize, 37), @as(usize, 40), @as(isize, -1), @as(isize, 1) },
    .{ .row_major, .conj_no_trans, @as(usize, 37), @as(usize, 37), @as(usize, 40), @as(isize, 2), @as(isize, -2) },
    .{ .col_major, .no_trans, @as(usize, 37), @as(usize, 37), @as(usize, 37), @as(isize, -2), @as(isize, -1) },
    .{ .col_major, .no_trans, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(isize, -1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 65), @as(usize, 65), @as(usize, 65), @as(isize, 1), @as(isize, -1) },

    // Thin matrices
    .{ .col_major, .no_trans, @as(usize, 10000), @as(usize, 4), @as(usize, 10000), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 3), @as(usize, 12000), @as(usize, 12000), @as(isize, 2), @as(isize, 1) },

    // Large arrays
    .{ .col_major, .no_trans, @as(usize, 2500), @as(usize, 2000), @as(usize, 2500), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 1000), @as(usize, 1500), @as(usize, 1500), @as(isize, 2), @as(isize, 1) },
    .{ .col_major, .conj_trans, @as(usize, 1500), @as(usize, 1000), @as(usize, 1500), @as(isize, -1), @as(isize, -1) },
    .{ .row_major, .conj_no_trans, @as(usize, 1500), @as(usize, 1000), @as(usize, 1000), @as(isize, 1), @as(isize, 2) },
    .{ .col_major, .trans, @as(usize, 1000), @as(usize, 1500), @as(usize, 1005), @as(isize, 3), @as(isize, 1) },
    .{ .col_major, .no_trans, @as(usize, 1500), @as(usize, 1000), @as(usize, 1500), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .conj_trans, @as(usize, 1000), @as(usize, 1500), @as(usize, 1500), @as(isize, -1), @as(isize, -1) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .no_trans, @as(usize, 1507), @as(usize, 1003), @as(usize, 1510), @as(isize, 1), @as(isize, 1) },
    .{ .row_major, .trans, @as(usize, 1003), @as(usize, 1507), @as(usize, 1510), @as(isize, 2), @as(isize, -1) },
};

test gemv {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeGemvTest(
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

fn executeGemvTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    transa: zsl.linalg.blas.Transpose,
    m: usize,
    n: usize,
    lda: usize,
    incx: isize,
    incy: isize,
    pool: *zsl.thread.Pool,
) !void {
    if (N == f64 and (transa == .conj_no_trans or transa == .conj_trans))
        return;

    const is_transposed = transa == .trans or transa == .conj_trans;

    const len_logical_x = if (is_transposed) m else n;
    const len_logical_y = if (is_transposed) n else m;

    const abs_incx = @abs(incx);
    const abs_incy = @abs(incy);

    const len_mem_x = if (len_logical_x == 0) 0 else 1 + (len_logical_x - 1) * abs_incx;
    const len_mem_y = if (len_logical_y == 0) 0 else 1 + (len_logical_y - 1) * abs_incy;

    const dim_a2 = if (layout == .col_major) n else m;
    const len_mem_a = lda * dim_a2;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const x = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x);

    const y_expected = try allocator.alloc(N, len_mem_y);
    defer allocator.free(y_expected);

    const y_actual = try allocator.alloc(N, len_mem_y);
    defer allocator.free(y_actual);
    const y_actual_parallel = try allocator.alloc(N, len_mem_y);
    defer allocator.free(y_actual_parallel);

    const alpha: N = tzsl.randomNumber(N, rand);
    const beta: N = tzsl.randomNumber(N, rand);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        y_actual[i] = elem.*;
        y_actual_parallel[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dgemv(
            layout.toInt(c_int),
            transa.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            beta,
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
        )
    else
        zsl.linalg.cblas.zgemv(
            layout.toInt(c_int),
            transa.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            &beta,
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
        );

    gemv(
        layout,
        transa,
        m,
        n,
        alpha,
        a.ptr,
        lda,
        x.ptr,
        incx,
        beta,
        y_actual.ptr,
        incy,
    ) catch |e| {
        std.debug.print("\n\tGEMV Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | transa: {s} | m: {} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), @tagName(layout), @tagName(transa), m, n, incx, incy });
        return e;
    };

    gemvParallel(
        layout,
        transa,
        m,
        n,
        alpha,
        a.ptr,
        lda,
        x.ptr,
        incx,
        beta,
        y_actual_parallel.ptr,
        incy,
        pool,
    ) catch |e| {
        std.debug.print("\n\tGEMV Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | transa: {s} | m: {} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), @tagName(layout), @tagName(transa), m, n, incx, incy });
        return e;
    };

    if (len_logical_y == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, len_logical_x);
    const abs_tol = 1e-10;
    for (y_expected, y_actual, y_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            const diff_parallel = zsl.float.abs(expected - actual_parallel);
            if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
                (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
            {
                std.debug.print("\n\tGEMV Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | transa: {s} | m: {} | n: {} | incx: {} | incy: {} | index: {}\n", .{ @tagName(layout), @tagName(transa), m, n, incx, incy, i });
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
                std.debug.print("\n\tGEMV Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | transa: {s} | m: {} | n: {} | incx: {} | incy: {} | index: {}\n", .{ @tagName(layout), @tagName(transa), m, n, incx, incy, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });
                std.debug.print("Actual (zsl, pool): {d} + {d}i\n", .{ actual_parallel.re, actual_parallel.im });
                std.debug.print("Diff (pool):        {d} + {d}i\n", .{ diff_parallel_re, diff_parallel_im });

                return error.TestFailed;
            }
        }
    }
}
