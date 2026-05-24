const std = @import("std");

const zsl = @import("zsl");
const nrm2 = zsl.linalg.blas.nrm2;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1), null, null },
    .{ @as(usize, 0), @as(isize, -1), null, null },

    // Single element
    .{ @as(usize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 1), @as(isize, 2), null, null },
    .{ @as(usize, 1), @as(isize, -1), null, null },

    // unroll aligned
    .{ @as(usize, 32), @as(isize, 1), null, null },
    .{ @as(usize, 64), @as(isize, 1), null, null },

    // unroll unaligned / remainder loops
    .{ @as(usize, 37), @as(isize, 1), null, null },
    .{ @as(usize, 69), @as(isize, 1), null, null },

    // Strided
    .{ @as(usize, 32), @as(isize, 2), null, null },
    .{ @as(usize, 65), @as(isize, 3), null, null },
    .{ @as(usize, 32), @as(isize, -1), null, null },
    .{ @as(usize, 65), @as(isize, -2), null, null },

    // Lowering parallel_threshold to test threading logic on small data chunks
    .{ @as(usize, 64), @as(isize, 1), null, @as(usize, 10) },
    .{ @as(usize, 65), @as(isize, 1), null, @as(usize, 8) },

    // Forced multithreaded
    .{ @as(usize, 65), @as(isize, 2), @as(usize, 2), null },
    .{ @as(usize, 65), @as(isize, -1), @as(usize, 4), null },
    .{ @as(usize, 65), @as(isize, -3), @as(usize, 2), null },

    // Default threading behavior
    .{ @as(usize, 1_500_000), @as(isize, 1), null, null },
    .{ @as(usize, 1_500_000), @as(isize, 2), null, null },
    .{ @as(usize, 1_500_000), @as(isize, -1), null, null },

    // Forced single-threaded fallback on large arrays
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(usize, 1), null },
    .{ @as(usize, 1_500_000), @as(isize, 3), @as(usize, 1), null },
    .{ @as(usize, 1_500_000), @as(isize, -2), @as(usize, 1), null },

    // Explicit high thread count
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(usize, 8), null },
    .{ @as(usize, 1_500_000), @as(isize, -1), @as(usize, 4), null },

    // Prime-like large numbers to ensure threads get unequal chunks
    .{ @as(usize, 1_500_007), @as(isize, 1), @as(usize, 4), null },
    .{ @as(usize, 1_500_007), @as(isize, 2), @as(usize, 0), null },
    .{ @as(usize, 1_500_007), @as(isize, -1), @as(usize, 2), null },
};

test nrm2 {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeNrm2Test(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                combo[3],
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
    num_threads: ?usize,
    parallel_threshold: ?usize,
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
        std.debug.print("\n\tNRM2 Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });

        return e;
    };

    if (n == 0 or expected == 0.0) {
        try std.testing.expectEqual(@as(f64, 0.0), actual);

        return;
    }

    const diff = zsl.float.abs(expected - actual);

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    if (diff > abs_tol and diff > expected * rel_tol) {
        std.debug.print("\n\tNRM2 Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });
        std.debug.print("Expected (CBLAS): {d}\n", .{expected});
        std.debug.print("Actual (zsl):     {d}\n", .{actual});
        std.debug.print("Diff:             {d}\n", .{diff});

        return error.TestFailed;
    }
}
