const std = @import("std");

const zsl = @import("zsl");
const ger = zsl.linalg.blas.ger;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, @as(usize, 0), @as(usize, 0), @as(isize, 1), @as(isize, 1), @as(usize, 1), null, null },
    .{ .row_major, @as(usize, 0), @as(usize, 10), @as(isize, -1), @as(isize, 2), @as(usize, 10), null, null },
    .{ .col_major, @as(usize, 10), @as(usize, 0), @as(isize, 1), @as(isize, -1), @as(usize, 10), null, null },

    // Single element
    .{ .row_major, @as(usize, 1), @as(usize, 1), @as(isize, 1), @as(isize, 1), @as(usize, 1), null, null },

    // Square aligned
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 32), null, null },
    .{ .row_major, @as(usize, 64), @as(usize, 64), @as(isize, 1), @as(isize, 1), @as(usize, 64), null, null },

    // Unaligned tall
    .{ .col_major, @as(usize, 69), @as(usize, 37), @as(isize, 1), @as(isize, 1), @as(usize, 69), null, null },

    // Unaligned wide
    .{ .row_major, @as(usize, 37), @as(usize, 69), @as(isize, 1), @as(isize, 1), @as(usize, 69), null, null },

    // Col major requires lda >= m
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 40), null, null },
    .{ .col_major, @as(usize, 37), @as(usize, 16), @as(isize, 1), @as(isize, 1), @as(usize, 45), null, null },

    // Row major requires lda >= n
    .{ .row_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 40), null, null },
    .{ .row_major, @as(usize, 16), @as(usize, 37), @as(isize, 1), @as(isize, 1), @as(usize, 45), null, null },

    // Strided
    .{ .col_major, @as(usize, 32), @as(usize, 32), @as(isize, 2), @as(isize, 1), @as(usize, 32), null, null },
    .{ .row_major, @as(usize, 32), @as(usize, 32), @as(isize, 1), @as(isize, 3), @as(usize, 32), null, null },
    .{ .col_major, @as(usize, 37), @as(usize, 37), @as(isize, -1), @as(isize, 1), @as(usize, 40), null, null },
    .{ .row_major, @as(usize, 37), @as(usize, 37), @as(isize, 2), @as(isize, -2), @as(usize, 40), null, null },
    .{ .col_major, @as(usize, 37), @as(usize, 37), @as(isize, -2), @as(isize, -1), @as(usize, 37), null, null },

    // Lowering parallel_threshold to test threading logic on small data chunks
    .{ .col_major, @as(usize, 64), @as(usize, 64), @as(isize, 1), @as(isize, 1), @as(usize, 64), null, @as(usize, 10) },
    .{ .row_major, @as(usize, 65), @as(usize, 65), @as(isize, 1), @as(isize, 1), @as(usize, 65), null, @as(usize, 8) },

    // Thin matrices
    .{ .col_major, @as(usize, 10000), @as(usize, 4), @as(isize, 1), @as(isize, 1), @as(usize, 10000), null, @as(usize, 10000) },
    .{ .row_major, @as(usize, 3), @as(usize, 12000), @as(isize, 2), @as(isize, 1), @as(usize, 12000), null, @as(usize, 9000) },

    // Forced multithreading
    .{ .col_major, @as(usize, 65), @as(usize, 65), @as(isize, 2), @as(isize, -1), @as(usize, 65), @as(usize, 2), null },
    .{ .row_major, @as(usize, 65), @as(usize, 65), @as(isize, -1), @as(isize, 2), @as(usize, 70), @as(usize, 4), null },
    .{ .col_major, @as(usize, 65), @as(usize, 65), @as(isize, -3), @as(isize, -3), @as(usize, 65), @as(usize, 2), null },

    // Large arrays
    .{ .col_major, @as(usize, 2500), @as(usize, 2000), @as(isize, 1), @as(isize, 1), @as(usize, 2500), null, null },
    .{ .row_major, @as(usize, 1000), @as(usize, 1500), @as(isize, 2), @as(isize, 1), @as(usize, 1500), null, null },

    // Forced single-threaded fallback
    .{ .row_major, @as(usize, 1500), @as(usize, 1000), @as(isize, 1), @as(isize, 2), @as(usize, 1500), @as(usize, 1), null },
    .{ .col_major, @as(usize, 1000), @as(usize, 1500), @as(isize, 3), @as(isize, 1), @as(usize, 1005), @as(usize, 1), null },

    // Explicit high thread count
    .{ .col_major, @as(usize, 1500), @as(usize, 1000), @as(isize, 1), @as(isize, 1), @as(usize, 1500), @as(usize, 8), null },
    .{ .row_major, @as(usize, 1000), @as(usize, 1500), @as(isize, -1), @as(isize, -1), @as(usize, 1500), @as(usize, 4), null },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, @as(usize, 1507), @as(usize, 1003), @as(isize, 1), @as(isize, 1), @as(usize, 1510), @as(usize, 4), null },
    .{ .row_major, @as(usize, 1003), @as(usize, 1507), @as(isize, 2), @as(isize, -1), @as(usize, 1510), @as(usize, 0), null },
};

test ger {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeGerTest(
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
            );
        }
    }
}

fn executeGerTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.Layout,
    m: usize,
    n: usize,
    incx: isize,
    incy: isize,
    lda: usize,
    num_threads: ?usize,
    parallel_threshold: ?usize,
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

    const alpha: N = tzsl.randomNumber(N, rand);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (a_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        a_actual[i] = elem.*;
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
        zsl.linalg.cblas.zgeru(
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

    ger(
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
        if (num_threads) |nt|
            if (parallel_threshold) |pt|
                .{
                    .num_threads = nt,
                    .parallel_threshold = pt,
                }
            else
                .{ .num_threads = nt }
        else if (parallel_threshold) |pt|
            .{
                .parallel_threshold = pt,
            }
        else
            .{},
    ) catch |e| {
        std.debug.print("\n\tGER Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {}\n", .{ @typeName(N), @tagName(layout), m, n, incx, incy, lda });
        return e;
    };

    if (m == 0 or n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, @max(m, n));
    const abs_tol = 1e-10;
    for (a_expected, a_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tGER Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {} | index: {}\n", .{ @tagName(layout), m, n, incx, incy, lda, i });
                std.debug.print("Expected (CBLAS): {d}\n", .{expected});
                std.debug.print("Actual (zsl):     {d}\n", .{actual});
                std.debug.print("Diff:             {d}\n", .{diff});

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
                std.debug.print("\n\tGER Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | m: {} | n: {} | incx: {} | incy: {} | lda: {} | index: {}\n", .{ @tagName(layout), m, n, incx, incy, lda, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:             {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
