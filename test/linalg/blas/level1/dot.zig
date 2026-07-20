const std = @import("std");

const zsl = @import("zsl");
const dot = zsl.linalg.blas.dot;
const dotParallel = zsl.linalg.blas.dotParallel;

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

test dot {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeDotTest(
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

fn executeDotTest(
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

    const x = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x);

    const y = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    const expected = if (N == f64)
        zsl.linalg.cblas.ddot(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y.ptr,
            zsl.numeric.cast(isize, incy),
        )
    else blk: {
        var result: zsl.linalg.blas.Dot([*]zsl.cf64, [*]zsl.cf64) = undefined;
        zsl.linalg.cblas.zdotu_sub(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y.ptr,
            zsl.numeric.cast(isize, incy),
            &result,
        );

        break :blk result;
    };

    const actual = dot(
        n,
        x.ptr,
        incx,
        y.ptr,
        incy,
    ) catch |e| {
        std.debug.print("\n\tDOT Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    const actual_parallel = dotParallel(
        n,
        x.ptr,
        incx,
        y.ptr,
        incy,
        pool,
    ) catch |e| {
        std.debug.print("\n\tDOT Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    if (n == 0) {
        if (N == f64) {
            try std.testing.expectEqual(@as(f64, 0.0), actual);
            try std.testing.expectEqual(@as(f64, 0.0), actual_parallel);
        } else {
            try std.testing.expectEqual(@as(f64, 0.0), actual.re);
            try std.testing.expectEqual(@as(f64, 0.0), actual.im);
            try std.testing.expectEqual(@as(f64, 0.0), actual_parallel.re);
            try std.testing.expectEqual(@as(f64, 0.0), actual_parallel.im);
        }

        return;
    }

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    if (N == f64) {
        const diff = zsl.float.abs(expected - actual);
        const diff_parallel = zsl.float.abs(expected - actual_parallel);
        if ((diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) or
            (diff_parallel > abs_tol and diff_parallel > zsl.float.abs(expected) * rel_tol))
        {
            std.debug.print("\n\tDOT Test Failed\n", .{});
            std.debug.print("Type: f64 | n: {} | incx: {} | incy: {}\n", .{ n, incx, incy });
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

        if ((diff_re > abs_tol and diff_re > zsl.float.abs(expected.re) * rel_tol) or
            (diff_im > abs_tol and diff_im > zsl.float.abs(expected.im) * rel_tol) or
            (diff_parallel_re > abs_tol and diff_parallel_re > zsl.float.abs(expected.re) * rel_tol) or
            (diff_parallel_im > abs_tol and diff_parallel_im > zsl.float.abs(expected.im) * rel_tol))
        {
            std.debug.print("\n\tDOT Test Failed\n", .{});
            std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {}\n", .{ n, incx, incy });
            std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
            std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
            std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });
            std.debug.print("Actual (zsl, pool): {d} + {d}i\n", .{ actual_parallel.re, actual_parallel.im });
            std.debug.print("Diff (pool):        {d} + {d}i\n", .{ diff_parallel_re, diff_parallel_im });

            return error.TestFailed;
        }
    }
}
