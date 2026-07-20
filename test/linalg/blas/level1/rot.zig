const std = @import("std");

const zsl = @import("zsl");
const rot = zsl.linalg.blas.rot;
const rotParallel = zsl.linalg.blas.rotParallel;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 0), @as(isize, -1), @as(isize, 2) },

    // Single element
    .{ @as(usize, 1), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 1), @as(isize, 2), @as(isize, -1) },
    .{ @as(usize, 1), @as(isize, -1), @as(isize, 2) },

    // unroll aligned
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 64), @as(isize, 1), @as(isize, 1) },

    // unroll unaligned / remainder loops
    .{ @as(usize, 37), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 69), @as(isize, 1), @as(isize, 1) },

    // Strided
    .{ @as(usize, 32), @as(isize, 2), @as(isize, 1) },
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 2) },
    .{ @as(usize, 65), @as(isize, 3), @as(isize, 3) },
    .{ @as(usize, 32), @as(isize, -1), @as(isize, 1) },
    .{ @as(usize, 65), @as(isize, 1), @as(isize, -2) },
    .{ @as(usize, 65), @as(isize, -2), @as(isize, -1) },

    // Small n relative to pool worker count: each worker gets a tiny chunk
    .{ @as(usize, 64), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 65), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 65), @as(isize, 2), @as(isize, -1) },
    .{ @as(usize, 65), @as(isize, -1), @as(isize, 2) },
    .{ @as(usize, 65), @as(isize, -3), @as(isize, -3) },

    // Large n, typical case
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 1_500_000), @as(isize, 2), @as(isize, 1) },
    .{ @as(usize, 1_500_000), @as(isize, -1), @as(isize, -1) },
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 2) },
    .{ @as(usize, 1_500_000), @as(isize, 3), @as(isize, 1) },
    .{ @as(usize, 1_500_000), @as(isize, -2), @as(isize, 2) },
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 1_500_000), @as(isize, -1), @as(isize, -1) },

    // Prime-like large numbers to ensure workers get unequal chunks
    .{ @as(usize, 1_500_007), @as(isize, 1), @as(isize, 1) },
    .{ @as(usize, 1_500_007), @as(isize, 2), @as(isize, -1) },
    .{ @as(usize, 1_500_007), @as(isize, -1), @as(isize, 1) },
};

test rot {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeRotTest(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                &pool,
            );
        }
    }
}

fn executeRotTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    incy: isize,
    pool: *zsl.thread.Pool,
) !void {
    const abs_incx = @abs(incx);
    const abs_incy = @abs(incy);

    const x_expected = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x_expected);

    const x_actual = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x_actual);
    const x_actual_parallel = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x_actual_parallel);

    const y_expected = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_expected);

    const y_actual = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_actual);
    const y_actual_parallel = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_actual_parallel);

    const c = tzsl.randomNumber(zsl.meta.Real(N), rand);
    const s = tzsl.randomNumber(zsl.meta.Real(N), rand);

    for (x_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);

        x_actual[i] = elem.*;
        x_actual_parallel[i] = elem.*;
    }

    for (y_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);

        y_actual[i] = elem.*;
        y_actual_parallel[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.drot(
            zsl.numeric.cast(isize, n),
            x_expected.ptr,
            zsl.numeric.cast(isize, incx),
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
            c,
            s,
        )
    else
        zsl.linalg.cblas.zdrot(
            zsl.numeric.cast(isize, n),
            x_expected.ptr,
            zsl.numeric.cast(isize, incx),
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
            c,
            s,
        );

    rot(
        n,
        x_actual.ptr,
        incx,
        y_actual.ptr,
        incy,
        c,
        s,
    ) catch |e| {
        std.debug.print("\n\tROT Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    rotParallel(
        n,
        x_actual_parallel.ptr,
        incx,
        y_actual_parallel.ptr,
        incy,
        c,
        s,
        pool,
    ) catch |e| {
        std.debug.print("\n\tROT Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    if (n == 0)
        return;

    const rel_tol = 1e-14;
    const abs_tol = 1e-10;
    for (x_expected, x_actual, x_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            const diff_parallel = zsl.float.abs(expected - actual_parallel);
            if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
                (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
            {
                std.debug.print("\n\tROT Test Failed\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (x, CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (x, zsl):       {d}\n", .{actual});
                std.debug.print("Diff (x):              {d}\n", .{diff});
                std.debug.print("Actual (x, zsl, pool): {d}\n", .{actual_parallel});
                std.debug.print("Diff (x, pool):        {d}\n", .{diff_parallel});

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
                std.debug.print("\n\tROT Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (x, CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (x, zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff (z):              {d} + {d}i\n", .{ diff_re, diff_im });
                std.debug.print("Actual (x, zsl, pool): {d} + {d}i\n", .{ actual_parallel.re, actual_parallel.im });
                std.debug.print("Diff (x, pool):        {d} + {d}i\n", .{ diff_parallel_re, diff_parallel_im });

                return error.TestFailed;
            }
        }
    }

    for (y_expected, y_actual, y_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            const diff_parallel = zsl.float.abs(expected - actual_parallel);
            if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
                (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
            {
                std.debug.print("\n\tROT Test Failed\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (y, CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (y, zsl):       {d}\n", .{actual});
                std.debug.print("Diff (y):              {d}\n", .{diff});
                std.debug.print("Actual (y, zsl, pool): {d}\n", .{actual_parallel});
                std.debug.print("Diff (y, pool):        {d}\n", .{diff_parallel});

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
                std.debug.print("\n\tROT Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (y, CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (y, zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff (z):              {d} + {d}i\n", .{ diff_re, diff_im });
                std.debug.print("Actual (y, zsl, pool): {d} + {d}i\n", .{ actual_parallel.re, actual_parallel.im });
                std.debug.print("Diff (y, pool):        {d} + {d}i\n", .{ diff_parallel_re, diff_parallel_im });

                return error.TestFailed;
            }
        }
    }
}
