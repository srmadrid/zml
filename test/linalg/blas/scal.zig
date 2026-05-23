const std = @import("std");

const zsl = @import("zsl");
const scal = zsl.linalg.blas.scal;

const tzsl = @import("../../zsl.zig");

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

test scal {
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
            try executeScalTest(T, allocator, rand, n, incx, num_threads, parallel_threshold);
        }
    }
}

fn executeScalTest(
    comptime T: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    num_threads: ?usize,
    parallel_threshold: ?usize,
) !void {
    const abs_incx = @abs(incx);

    const len_x = if (n == 0) 0 else 1 + (n - 1) * abs_incx;

    const x_expected = try allocator.alloc(T, len_x);
    defer allocator.free(x_expected);

    const x_actual = try allocator.alloc(T, len_x);
    defer allocator.free(x_actual);

    var alpha: T = undefined;
    if (T == f64) {
        alpha = rand.float(f64) * 2.0 - 1.0;
    } else {
        alpha = .{
            .re = rand.float(f64) * 2.0 - 1.0,
            .im = rand.float(f64) * 2.0 - 1.0,
        };
    }

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

    if (T == f64) {
        zsl.linalg.cblas.dscal(zsl.numeric.cast(isize, n), alpha, @ptrCast(x_expected.ptr), zsl.numeric.cast(isize, abs_incx));
    } else {
        zsl.linalg.cblas.zscal(zsl.numeric.cast(isize, n), &alpha, @ptrCast(x_expected.ptr), zsl.numeric.cast(isize, abs_incx));
    }

    scal(
        n,
        alpha,
        x_actual.ptr,
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
        std.debug.print("\n\tSCAL Test Failed (Exception)\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(T), n, incx });
        return e;
    };

    if (n == 0) return;

    const rel_tol = 1e-12;
    const abs_tol = 1e-10;

    for (x_expected, x_actual, 0..) |expected, actual, i| {
        if (T == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tSCAL Test Failed (Data Mismatch)\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | index: {}\n", .{ n, incx, i });
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
                std.debug.print("\n\tSCAL Test Failed (Data Mismatch)\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | index: {}\n", .{ n, incx, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:             {d} + {d}i\n", .{ diff_re, diff_im });
                return error.TestFailed;
            }
        }
    }
}
