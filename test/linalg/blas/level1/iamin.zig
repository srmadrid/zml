const std = @import("std");

const zsl = @import("zsl");
const iamin = zsl.linalg.blas.iamin;
const iaminParallel = zsl.linalg.blas.iaminParallel;

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

test iamin {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeIaminTest(
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

fn executeIaminTest(
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

    if (n > 0) {
        const target_logical_idx = rand.intRangeLessThan(usize, 0, n);
        const mem_idx = if (incx < 0) (n - 1 - target_logical_idx) * abs_incx else target_logical_idx * abs_incx;

        if (N == f64) {
            x[mem_idx] = 0.0;
        } else {
            x[mem_idx] = .{ .re = 0.0, .im = 0.0 };
        }
    }

    const cblas_idx_raw = if (N == f64)
        zsl.linalg.cblas.idamin(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, abs_incx),
        )
    else
        zsl.linalg.cblas.izamin(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, abs_incx),
        );

    const expected = if (n == 0) @as(usize, 0) else if (incx < 0) n - 1 - zsl.numeric.cast(usize, cblas_idx_raw) else zsl.numeric.cast(usize, cblas_idx_raw);

    const actual = iamin(
        n,
        x.ptr,
        incx,
    ) catch |e| {
        std.debug.print("\n\tIAMIN Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });

        return e;
    };

    const actual_parallel = iaminParallel(
        n,
        x.ptr,
        incx,
        pool,
    ) catch |e| {
        std.debug.print("\n\tIAMIN Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });

        return e;
    };

    if (n == 0) {
        try std.testing.expectEqual(@as(usize, 0), actual);
        try std.testing.expectEqual(@as(usize, 0), actual_parallel);

        return;
    }

    if (expected != actual or expected != actual_parallel) {
        std.debug.print("\n\tIAMIN Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });
        std.debug.print("Expected (CBLAS mapped): {}\n", .{expected});
        std.debug.print("Actual (zsl):            {}\n", .{actual});
        std.debug.print("Actual (zsl, pool):      {}\n", .{actual_parallel});

        return error.TestFailed;
    }
}
