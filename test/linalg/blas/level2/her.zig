const std = @import("std");

const zsl = @import("zsl");
const her = zsl.linalg.blas.her;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, .upper, @as(usize, 0), @as(isize, 1), @as(usize, 1) },
    .{ .row_major, .lower, @as(usize, 0), @as(isize, -1), @as(usize, 10) },

    // Single element
    .{ .col_major, .upper, @as(usize, 1), @as(isize, 1), @as(usize, 1) },
    .{ .row_major, .lower, @as(usize, 1), @as(isize, -1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .upper, @as(usize, 32), @as(isize, 1), @as(usize, 32) },
    .{ .row_major, .lower, @as(usize, 64), @as(isize, 1), @as(usize, 64) },

    // Requires lda >= n
    .{ .col_major, .upper, @as(usize, 32), @as(isize, 1), @as(usize, 40) },
    .{ .row_major, .lower, @as(usize, 37), @as(isize, 1), @as(usize, 45) },

    // Strided
    .{ .col_major, .lower, @as(usize, 32), @as(isize, 2), @as(usize, 32) },
    .{ .row_major, .upper, @as(usize, 32), @as(isize, 3), @as(usize, 32) },
    .{ .col_major, .upper, @as(usize, 37), @as(isize, -1), @as(usize, 40) },
    .{ .row_major, .lower, @as(usize, 37), @as(isize, -2), @as(usize, 40) },
    .{ .col_major, .lower, @as(usize, 37), @as(isize, -2), @as(usize, 37) },
    .{ .col_major, .upper, @as(usize, 64), @as(isize, -1), @as(usize, 64) },
    .{ .row_major, .lower, @as(usize, 65), @as(isize, -7), @as(usize, 65) },
    .{ .col_major, .lower, @as(usize, 65), @as(isize, 2), @as(usize, 65) },
    .{ .row_major, .upper, @as(usize, 65), @as(isize, -1), @as(usize, 70) },
    .{ .col_major, .upper, @as(usize, 65), @as(isize, -3), @as(usize, 65) },

    // Large arrays
    .{ .col_major, .upper, @as(usize, 2000), @as(isize, 1), @as(usize, 2000) },
    .{ .row_major, .lower, @as(usize, 1500), @as(isize, 2), @as(usize, 1500) },
    .{ .row_major, .lower, @as(usize, 1500), @as(isize, 1), @as(usize, 1500) },
    .{ .col_major, .upper, @as(usize, 1000), @as(isize, 3), @as(usize, 1005) },
    .{ .col_major, .lower, @as(usize, 1500), @as(isize, 1), @as(usize, 1500) },
    .{ .row_major, .upper, @as(usize, 1000), @as(isize, -1), @as(usize, 1000) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .upper, @as(usize, 1507), @as(isize, 1), @as(usize, 1510) },
    .{ .row_major, .lower, @as(usize, 1003), @as(isize, 2), @as(usize, 1010) },
};

test her {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        try executeHerTest(
            zsl.cf64,
            allocator,
            rand,
            combo[0],
            combo[1],
            combo[2],
            combo[3],
            combo[4],
            &pool,
        );
    }
}

fn executeHerTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    uplo: zsl.matrix.Uplo,
    n: usize,
    incx: isize,
    lda: usize,
    pool: *zsl.thread.Pool,
) !void {
    _ = pool;

    const abs_incx = @abs(incx);

    const len_mem_x = if (n == 0) 0 else 1 + (n - 1) * abs_incx;
    const len_mem_a = lda * n;

    const x = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x);

    const a_expected = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_expected);

    const a_actual = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_actual);

    const alpha: f64 = tzsl.randomNumber(f64, rand);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (a_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        a_actual[i] = elem.*;
    }

    var diag_idx: usize = 0;
    while (diag_idx < n) : (diag_idx += 1) {
        a_expected[diag_idx * lda + diag_idx].im = 0;
        a_actual[diag_idx * lda + diag_idx].im = 0;
    }

    zsl.linalg.cblas.zher(
        layout.toInt(c_int),
        uplo.toInt(c_int),
        zsl.numeric.cast(isize, n),
        alpha,
        x.ptr,
        zsl.numeric.cast(isize, incx),
        a_expected.ptr,
        zsl.numeric.cast(isize, lda),
    );

    her(
        layout,
        uplo,
        n,
        alpha,
        x.ptr,
        incx,
        a_actual.ptr,
        lda,
    ) catch |e| {
        std.debug.print("\n\tHER Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | n: {} | incx: {} | lda: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), n, incx, lda });
        return e;
    };

    if (n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    for (a_expected, a_actual, 0..) |expected, actual, i| {
        const diff_re = zsl.float.abs(expected.re - actual.re);
        const diff_im = zsl.float.abs(expected.im - actual.im);

        const mag_re = zsl.float.abs(expected.re);
        const mag_im = zsl.float.abs(expected.im);

        if ((diff_re > abs_tol and diff_re > mag_re * rel_tol) or
            (diff_im > abs_tol and diff_im > mag_im * rel_tol))
        {
            std.debug.print("\n\tHER Test Failed\n", .{});
            std.debug.print("Type: zsl.cf64 | layout: {s} | uplo: {s} | n: {} | incx: {} | lda: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), n, incx, lda, i });
            std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
            std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
            std.debug.print("Diff:             {d} + {d}i\n", .{ diff_re, diff_im });

            return error.TestFailed;
        }
    }
}
