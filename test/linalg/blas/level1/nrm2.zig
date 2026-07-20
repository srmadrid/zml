const std = @import("std");

const zsl = @import("zsl");
const nrm2 = zsl.linalg.blas.nrm2;
const nrm2Parallel = zsl.linalg.blas.nrm2Parallel;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1) },
    .{ @as(usize, 0), @as(isize, -1) },

    // Single element
    .{ @as(usize, 1), @as(isize, 1) },
    .{ @as(usize, 1), @as(isize, 2) },
    .{ @as(usize, 1), @as(isize, -1) },

    // Unroll aligned
    .{ @as(usize, 32), @as(isize, 1) },
    .{ @as(usize, 64), @as(isize, 1) },

    // Unroll unaligned / remainder loops
    .{ @as(usize, 37), @as(isize, 1) },
    .{ @as(usize, 69), @as(isize, 1) },

    // Strided
    .{ @as(usize, 32), @as(isize, 2) },
    .{ @as(usize, 65), @as(isize, 3) },
    .{ @as(usize, 32), @as(isize, -1) },
    .{ @as(usize, 65), @as(isize, -2) },

    // Small n relative to pool worker count: each worker gets a tiny chunk
    .{ @as(usize, 65), @as(isize, 1) },
    .{ @as(usize, 65), @as(isize, 2) },
    .{ @as(usize, 65), @as(isize, -1) },
    .{ @as(usize, 65), @as(isize, -3) },

    // Large n, typical case
    .{ @as(usize, 1_500_000), @as(isize, 1) },
    .{ @as(usize, 1_500_000), @as(isize, 2) },
    .{ @as(usize, 1_500_000), @as(isize, -1) },
    .{ @as(usize, 1_500_000), @as(isize, 3) },
    .{ @as(usize, 1_500_000), @as(isize, -2) },

    // Prime-like large numbers to ensure workers get unequal chunks
    .{ @as(usize, 1_500_007), @as(isize, 1) },
    .{ @as(usize, 1_500_007), @as(isize, 2) },
    .{ @as(usize, 1_500_007), @as(isize, -1) },
};

test nrm2 {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeNrm2Test(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                &pool,
            );
        }
    }
}

fn executeNrm2Test(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    pool: *zsl.thread.Pool,
) !void {
    const abs_incx = @abs(incx);

    const x = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    const expected: f64 = if (N == f64)
        zsl.linalg.cblas.dnrm2(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, abs_incx),
        )
    else
        zsl.linalg.cblas.dznrm2(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, abs_incx),
        );

    const actual = nrm2(
        n,
        x.ptr,
        incx,
    ) catch |e| {
        std.debug.print("\n\tNRM2 Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });

        return e;
    };

    const actual_parallel = nrm2Parallel(
        n,
        x.ptr,
        incx,
        pool,
    ) catch |e| {
        std.debug.print("\n\tNRM2 Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });

        return e;
    };

    if (n == 0 or expected == 0.0) {
        try std.testing.expectEqual(@as(f64, 0.0), actual);
        try std.testing.expectEqual(@as(f64, 0.0), actual_parallel);

        return;
    }

    const diff = zsl.float.abs(expected - actual);
    const diff_parallel = zsl.float.abs(expected - actual_parallel);

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    if ((diff > abs_tol and diff > expected * rel_tol) or
        (diff_parallel > abs_tol and diff_parallel > expected * rel_tol))
    {
        std.debug.print("\n\tNRM2 Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });
        std.debug.print("Expected (CBLAS):   {d}\n", .{expected});
        std.debug.print("Actual (zsl):       {d}\n", .{actual});
        std.debug.print("Diff:               {d}\n", .{diff});
        std.debug.print("Actual (zsl, pool): {d}\n", .{actual_parallel});
        std.debug.print("Diff (pool):        {d}\n", .{diff_parallel});

        return error.TestFailed;
    }
}
