const std = @import("std");

const zsl = @import("zsl");
const dot = zsl.linalg.blas.dot;

const tzsl = @import("../../zsl.zig");

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

test dot {
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
            try executeDotTest(T, allocator, rand, n, incx, incy, num_threads, parallel_threshold);
        }
    }
}

fn executeDotTest(
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

    const x = try allocator.alloc(T, len_x);
    defer allocator.free(x);

    const y = try allocator.alloc(T, len_y);
    defer allocator.free(y);

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

    for (y) |*elem| {
        if (T == f64) {
            elem.* = rand.float(f64) * 20.0 - 10.0;
        } else {
            elem.* = .{
                .re = rand.float(f64) * 20.0 - 10.0,
                .im = rand.float(f64) * 20.0 - 10.0,
            };
        }
    }

    const expected = if (T == f64)
        zsl.linalg.cblas.ddot(zsl.numeric.cast(isize, n), @ptrCast(x.ptr), zsl.numeric.cast(isize, incx), @ptrCast(y.ptr), zsl.numeric.cast(isize, incy))
    else blk: {
        var result: zsl.linalg.blas.Dot([*]zsl.cf64, [*]zsl.cf64) = undefined;
        zsl.linalg.cblas.zdotu_sub(zsl.numeric.cast(isize, n), @ptrCast(x.ptr), zsl.numeric.cast(isize, incx), @ptrCast(y.ptr), zsl.numeric.cast(isize, incy), &result);
        break :blk result;
    };

    const actual = dot(
        n,
        x.ptr,
        incx,
        y.ptr,
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
        std.debug.print("\n\tDOT Test Failed (Exception)\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(T), n, incx, incy });
        return e;
    };

    if (n == 0) {
        if (T == f64) {
            try std.testing.expectEqual(@as(f64, 0.0), actual);
        } else {
            try std.testing.expectEqual(@as(f64, 0.0), actual.re);
            try std.testing.expectEqual(@as(f64, 0.0), actual.im);
        }
        return;
    }

    const rel_tol = 1e-13 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    if (T == f64) {
        const diff = zsl.float.abs(expected - actual);
        if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
            std.debug.print("\n\tDOT Test Failed\n", .{});
            std.debug.print("Type: f64 | n: {} | incx: {} | incy: {}\n", .{ n, incx, incy });
            std.debug.print("Expected (CBLAS): {d}\n", .{expected});
            std.debug.print("Actual (zsl):     {d}\n", .{actual});
            std.debug.print("Diff:             {d}\n", .{diff});
            return error.TestFailed;
        }
    } else {
        const diff_re = zsl.float.abs(expected.re - actual.re);
        const diff_im = zsl.float.abs(expected.im - actual.im);

        if ((diff_re > abs_tol and diff_re > zsl.float.abs(expected.re) * rel_tol) or
            (diff_im > abs_tol and diff_im > zsl.float.abs(expected.im) * rel_tol))
        {
            std.debug.print("\n\tDOT Test Failed\n", .{});
            std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {}\n", .{ n, incx, incy });
            std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
            std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
            std.debug.print("Diff:             {d} + {d}i\n", .{ diff_re, diff_im });
            return error.TestFailed;
        }
    }
}
