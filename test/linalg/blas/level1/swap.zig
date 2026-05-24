const std = @import("std");

const zsl = @import("zsl");
const swap = zsl.linalg.blas.swap;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 0), @as(isize, -1), @as(isize, 2), null, null },

    // Single element
    .{ @as(usize, 1), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 1), @as(isize, 2), @as(isize, -1), null, null },
    .{ @as(usize, 1), @as(isize, -1), @as(isize, 2), null, null },

    // SIMD aligned
    .{ @as(usize, 16), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 1), null, null },

    // SIMD unaligned / remainder loops
    .{ @as(usize, 7), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 33), @as(isize, 1), @as(isize, 1), null, null },

    // Mixed Strides
    .{ @as(usize, 16), @as(isize, 2), @as(isize, 1), null, null },
    .{ @as(usize, 16), @as(isize, 1), @as(isize, 2), null, null },
    .{ @as(usize, 33), @as(isize, 3), @as(isize, 3), null, null },

    // Reverse strides
    .{ @as(usize, 16), @as(isize, -1), @as(isize, 1), null, null },
    .{ @as(usize, 33), @as(isize, 1), @as(isize, -2), null, null },
    .{ @as(usize, 33), @as(isize, -2), @as(isize, -1), null, null },

    // Lowering parallel_threshold to test threading logic on small data chunks
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 2), @as(usize, 10) },
    .{ @as(usize, 33), @as(isize, 1), @as(isize, 1), @as(usize, 4), @as(usize, 8) },

    // Forced multithreaded + strided
    .{ @as(usize, 33), @as(isize, 2), @as(isize, -1), @as(usize, 2), @as(usize, 10) },
    .{ @as(usize, 33), @as(isize, -1), @as(isize, 2), @as(usize, 4), @as(usize, 8) },
    .{ @as(usize, 33), @as(isize, -3), @as(isize, -3), @as(usize, 2), @as(usize, 10) },

    // Default threading behavior (opts.num_threads = 0)
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 1_500_000), @as(isize, 2), @as(isize, 1), null, null },
    .{ @as(usize, 1_500_000), @as(isize, -1), @as(isize, -1), null, null },

    // Forced single-threaded fallback on large arrays
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 2), @as(usize, 1), null },
    .{ @as(usize, 1_500_000), @as(isize, 3), @as(isize, 1), @as(usize, 1), null },
    .{ @as(usize, 1_500_000), @as(isize, -2), @as(isize, 2), @as(usize, 1), null },

    // Explicit high thread count
    .{ @as(usize, 1_500_000), @as(isize, 1), @as(isize, 1), @as(usize, 8), null },
    .{ @as(usize, 1_500_000), @as(isize, -1), @as(isize, -1), @as(usize, 4), null },

    // Prime-like large numbers to ensure threads get unequal chunks
    .{ @as(usize, 1_500_007), @as(isize, 1), @as(isize, 1), @as(usize, 4), null },
    .{ @as(usize, 1_500_007), @as(isize, 2), @as(isize, -1), @as(usize, 0), null },
    .{ @as(usize, 1_500_007), @as(isize, -1), @as(isize, 1), @as(usize, 2), null },
};

test swap {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |combo| {
        const n = combo[0];
        const incx = combo[1];
        const incy = combo[2];
        const num_threads = combo[3];
        const parallel_threshold = combo[4];

        inline for (.{ f64, zsl.cf64 }) |T| {
            try executeSwapTest(T, allocator, rand, n, incx, incy, num_threads, parallel_threshold);
        }
    }
}

fn executeSwapTest(
    comptime T: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    incy: isize,
    num_threads: ?usize,
    parallel_threshold: ?usize,
) !void {
    const abs_incx = @abs(incx);
    const abs_incy = @abs(incy);

    const len_x = if (n == 0) 0 else 1 + (n - 1) * abs_incx;
    const len_y = if (n == 0) 0 else 1 + (n - 1) * abs_incy;

    const x_expected = try allocator.alloc(T, len_x);
    defer allocator.free(x_expected);

    const x_actual = try allocator.alloc(T, len_x);
    defer allocator.free(x_actual);

    const y_expected = try allocator.alloc(T, len_y);
    defer allocator.free(y_expected);

    const y_actual = try allocator.alloc(T, len_y);
    defer allocator.free(y_actual);

    for (x_expected, 0..) |*elem, i| {
        if (T == f64) {
            elem.* = rand.float(f64) * 20.0 - 10.0;
        } else {
            elem.* = .{
                .re = rand.float(f64) * 20.0 - 10.0,
                .im = rand.float(f64) * 20.0 - 10.0,
            };
        }

        x_actual[i] = elem.*;
    }

    for (y_expected, 0..) |*elem, i| {
        if (T == f64) {
            elem.* = rand.float(f64) * 20.0 - 10.0;
        } else {
            elem.* = .{
                .re = rand.float(f64) * 20.0 - 10.0,
                .im = rand.float(f64) * 20.0 - 10.0,
            };
        }

        y_actual[i] = elem.*;
    }

    if (T == f64) {
        zsl.linalg.cblas.dswap(zsl.numeric.cast(isize, n), @ptrCast(x_expected.ptr), zsl.numeric.cast(isize, incx), @ptrCast(y_expected.ptr), zsl.numeric.cast(isize, incy));
    } else {
        zsl.linalg.cblas.zswap(zsl.numeric.cast(isize, n), @ptrCast(x_expected.ptr), zsl.numeric.cast(isize, incx), @ptrCast(y_expected.ptr), zsl.numeric.cast(isize, incy));
    }

    swap(
        n,
        x_actual.ptr,
        incx,
        y_actual.ptr,
        incy,
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
        std.debug.print("\n\tSWAP Test Failed (Exception)\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(T), n, incx, incy });
        return e;
    };

    if (n == 0)
        return;

    for (x_expected, x_actual, 0..) |expected, actual, i| {
        if (T == f64) {
            if (expected != actual) {
                std.debug.print("\n\tSWAP Test Failed (X Data Mismatch)\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d}\n", .{expected});
                std.debug.print("Actual (zsl):     {d}\n", .{actual});
                return error.TestFailed;
            }
        } else {
            if (expected.re != actual.re or expected.im != actual.im) {
                std.debug.print("\n\tSWAP Test Failed (X Data Mismatch)\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
                return error.TestFailed;
            }
        }
    }

    for (y_expected, y_actual, 0..) |expected, actual, i| {
        if (T == f64) {
            if (expected != actual) {
                std.debug.print("\n\tSWAP Test Failed (Y Data Mismatch)\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d}\n", .{expected});
                std.debug.print("Actual (zsl):     {d}\n", .{actual});
                return error.TestFailed;
            }
        } else {
            if (expected.re != actual.re or expected.im != actual.im) {
                std.debug.print("\n\tSWAP Test Failed (Y Data Mismatch)\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
                return error.TestFailed;
            }
        }
    }
}
