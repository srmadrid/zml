const std = @import("std");

const zsl = @import("zsl");
const gerc = zsl.linalg.blas.gerc;
const gercParallel = zsl.linalg.blas.gercParallel;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, @as(usize, 0), @as(usize, 0), @as(isize, 1), @as(isize, 1), @as(usize, 1) },
    .{ .row_major, @as(usize, 0), @as(usize, 10), @as(isize, -1), @as(isize, 2), @as(usize, 10) },
    .{ .col_major, @as(usize, 10), @as(usize, 0), @as(isize, 1), @as(isize, -1), @as(usize, 10) },

    // Single element
    .{ .row_major, @as(usize, 1), @as(usize, 1), @as(isize, 1), @as(isize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 32) },
    .{ .row_major, @as(usize, 64), @as(usize, 64), @as(isize, 1), @as(isize, 1), @as(usize, 64) },

    // Unaligned tall
    .{ .col_major, @as(usize, 69), @as(usize, 37), @as(isize, 1), @as(isize, 1), @as(usize, 69) },

    // Unaligned wide
    .{ .row_major, @as(usize, 37), @as(usize, 69), @as(isize, 1), @as(isize, 1), @as(usize, 69) },

    // Col major requires lda >= m
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 40) },
    .{ .col_major, @as(usize, 37), @as(usize, 16), @as(isize, 1), @as(isize, 1), @as(usize, 45) },

    // Row major requires lda >= n
    .{ .row_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 40) },
    .{ .row_major, @as(usize, 16), @as(usize, 37), @as(isize, 1), @as(isize, 1), @as(usize, 45) },

    // Strided
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 2), @as(isize, 1), @as(usize, 32) },
    .{ .row_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 3), @as(usize, 32) },
    .{ .col_major, @as(usize, 37), @as(usize, 37), @as(isize, -1), @as(isize, 1), @as(usize, 40) },
    .{ .row_major, @as(usize, 37), @as(usize, 37), @as(isize, 2), @as(isize, -2), @as(usize, 40) },
    .{ .col_major, @as(usize, 37), @as(usize, 37), @as(isize, -2), @as(isize, -1), @as(usize, 37) },
    .{ .col_major, @as(usize, 64), @as(usize, 64), @as(isize, -1), @as(isize, 1), @as(usize, 64) },
    .{ .row_major, @as(usize, 65), @as(usize, 65), @as(isize, 1), @as(isize, -1), @as(usize, 65) },

    // Thin matrices
    .{ .col_major, @as(usize, 10000), @as(usize, 4), @as(isize, 1), @as(isize, 1), @as(usize, 10000) },
    .{ .row_major, @as(usize, 3), @as(usize, 12000), @as(isize, 2), @as(isize, 1), @as(usize, 12000) },

    // Large arrays
    .{ .col_major, @as(usize, 2500), @as(usize, 2000), @as(isize, 1), @as(isize, 1), @as(usize, 2500) },
    .{ .row_major, @as(usize, 1000), @as(usize, 1500), @as(isize, 2), @as(isize, 1), @as(usize, 1500) },
    .{ .col_major, @as(usize, 1500), @as(usize, 1000), @as(isize, -1), @as(isize, -1), @as(usize, 1500) },
    .{ .row_major, @as(usize, 1500), @as(usize, 1000), @as(isize, 1), @as(isize, 2), @as(usize, 1000) },
    .{ .col_major, @as(usize, 1000), @as(usize, 1500), @as(isize, 3), @as(isize, 1), @as(usize, 1005) },
    .{ .col_major, @as(usize, 1500), @as(usize, 1000), @as(isize, 1), @as(isize, 1), @as(usize, 1500) },
    .{ .row_major, @as(usize, 1000), @as(usize, 1500), @as(isize, -1), @as(isize, -1), @as(usize, 1500) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, @as(usize, 1507), @as(usize, 1003), @as(isize, 1), @as(isize, 1), @as(usize, 1510) },
    .{ .row_major, @as(usize, 1003), @as(usize, 1507), @as(isize, 2), @as(isize, -1), @as(usize, 1510) },
};

test gerc {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeGercTest(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                combo[3],
                combo[4],
                combo[5],
                &pool,
            );
        }
    }
}

fn executeGercTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    m: usize,
    n: usize,
    incx: isize,
    incy: isize,
    lda: usize,
    pool: *zsl.thread.Pool,
) !void {
    const abs_incx = @abs(incx);
    const abs_incy = @abs(incy);

    const len_mem_x = if (m == 0) 0 else 1 + (m - 1) * abs_incx;
    const len_mem_y = if (n == 0) 0 else 1 + (n - 1) * abs_incy;

    const dim_a2 = if (layout == .col_major) n else m;
    const len_mem_a = lda * dim_a2;

    const x = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x);

    const y = try allocator.alloc(N, len_mem_y);
    defer allocator.free(y);

    const a_expected = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_expected);

    const a_actual = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_actual);
    const a_actual_parallel = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_actual_parallel);

    const alpha: N = tzsl.randomNumber(N, rand);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (a_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        a_actual[i] = elem.*;
        a_actual_parallel[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dger(
            layout.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            alpha,
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y.ptr,
            zsl.numeric.cast(isize, incy),
            a_expected.ptr,
            zsl.numeric.cast(isize, lda),
        )
    else
        zsl.linalg.cblas.zgerc(
            layout.toInt(c_int),
            zsl.numeric.cast(isize, m),
            zsl.numeric.cast(isize, n),
            &alpha,
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y.ptr,
            zsl.numeric.cast(isize, incy),
            a_expected.ptr,
            zsl.numeric.cast(isize, lda),
        );

    gerc(
        layout,
        m,
        n,
        alpha,
        x.ptr,
        incx,
        y.ptr,
        incy,
        a_actual.ptr,
        lda,
    ) catch |e| {
        std.debug.print("\n\tGERC Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {}\n", .{ @typeName(N), @tagName(layout), m, n, incx, incy, lda });
        return e;
    };

    gercParallel(
        layout,
        m,
        n,
        alpha,
        x.ptr,
        incx,
        y.ptr,
        incy,
        a_actual_parallel.ptr,
        lda,
        pool,
    ) catch |e| {
        std.debug.print("\n\tGERC Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {}\n", .{ @typeName(N), @tagName(layout), m, n, incx, incy, lda });
        return e;
    };

    if (m == 0 or n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, @max(m, n));
    const abs_tol = 1e-10;
    for (a_expected, a_actual, a_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            const diff_parallel = zsl.float.abs(expected - actual_parallel);
            if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
                (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
            {
                std.debug.print("\n\tGERC Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {} | index: {}\n", .{ @tagName(layout), m, n, incx, incy, lda, i });
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
                std.debug.print("\n\tGERC Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {} | index: {}\n", .{ @tagName(layout), m, n, incx, incy, lda, i });
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
