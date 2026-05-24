const std = @import("std");

const zsl = @import("zsl");
const iamax = zsl.linalg.blas.iamax;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1), null, null },
    .{ @as(usize, 0), @as(isize, -1), null, null },

    // Single element
    .{ @as(usize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 1), @as(isize, 2), null, null },
    .{ @as(usize, 1), @as(isize, -1), null, null },

    // SIMD aligned
    .{ @as(usize, 16), @as(isize, 1), null, null },
    .{ @as(usize, 32), @as(isize, 1), null, null },

    // SIMD unaligned / remainder loops
    .{ @as(usize, 7), @as(isize, 1), null, null },
    .{ @as(usize, 33), @as(isize, 1), null, null },

    // Strided dense (incx > 1)
    .{ @as(usize, 16), @as(isize, 2), null, null },
    .{ @as(usize, 33), @as(isize, 3), null, null },

    // Reverse strides (incx < 0)
    .{ @as(usize, 16), @as(isize, -1), null, null },
    .{ @as(usize, 33), @as(isize, -2), null, null },

    // Lowering parallel_threshold to test threading logic on small data chunks
    .{ @as(usize, 32), @as(isize, 1), @as(usize, 2), @as(usize, 10) },
    .{ @as(usize, 33), @as(isize, 1), @as(usize, 4), @as(usize, 8) },

    // Forced multithreaded + strided
    .{ @as(usize, 33), @as(isize, 2), @as(usize, 2), @as(usize, 10) },
    .{ @as(usize, 33), @as(isize, -1), @as(usize, 4), @as(usize, 8) },
    .{ @as(usize, 33), @as(isize, -3), @as(usize, 2), @as(usize, 10) },

    // Default threading behavior (opts.num_threads = 0)
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

test iamax {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |combo| {
        const n = combo[0];
        const incx = combo[1];
        const num_threads = combo[2];
        const parallel_threshold = combo[3];

        inline for (.{ f64, zsl.cf64 }) |T| {
            try executeIamaxTest(T, allocator, rand, n, incx, num_threads, parallel_threshold);
        }
    }
}

fn executeIamaxTest(
    comptime T: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    num_threads: ?usize,
    parallel_threshold: ?usize,
) !void {
    const abs_incx = @abs(incx);

    const x = try allocator.alloc(T, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x);

    for (x) |*elem| {
        if (T == f64) {
            elem.* = rand.float(f64) * 20.0 - 10.0;
        } else {
            elem.* = .{
                .re = rand.float(f64) * 20.0 - 10.0,
                .im = rand.float(f64) * 20.0 - 10.0,
            };
        }
    }

    if (n > 0) {
        const target_logical_idx = rand.intRangeLessThan(usize, 0, n);
        const mem_idx = if (incx < 0) (n - 1 - target_logical_idx) * abs_incx else target_logical_idx * abs_incx;

        if (T == f64) {
            x[mem_idx] = 10000.0;
        } else {
            x[mem_idx] = .{ .re = 10000.0, .im = 10000.0 };
        }
    }

    const cblas_idx_raw = if (T == f64)
        zsl.linalg.cblas.idamax(zsl.numeric.cast(isize, n), @ptrCast(x.ptr), zsl.numeric.cast(isize, abs_incx))
    else
        zsl.linalg.cblas.izamax(zsl.numeric.cast(isize, n), @ptrCast(x.ptr), zsl.numeric.cast(isize, abs_incx));

    const expected = if (n == 0) @as(usize, 0) else if (incx < 0) n - 1 - zsl.numeric.cast(usize, cblas_idx_raw) else zsl.numeric.cast(usize, cblas_idx_raw);

    const actual = iamax(
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
        std.debug.print("\n\tIAMAX Test Failed (Exception)\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(T), n, incx });

        return e;
    };

    if (n == 0) {
        try std.testing.expectEqual(@as(usize, 0), actual);
        return;
    }

    if (expected != actual) {
        std.debug.print("\n\tIAMAX Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(T), n, incx });
        std.debug.print("Expected (CBLAS mapped): {}\n", .{expected});
        std.debug.print("Actual (zsl):            {}\n", .{actual});

        return error.TestFailed;
    }
}
