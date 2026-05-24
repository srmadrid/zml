const std = @import("std");

const zsl = @import("zsl");
const scal = zsl.linalg.blas.scal;

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
    .{ @as(usize, 16), @as(isize, 1), null, null },
    .{ @as(usize, 32), @as(isize, 1), null, null },

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
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeScalTest(
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

fn executeScalTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    n: usize,
    incx: isize,
    num_threads: ?usize,
    parallel_threshold: ?usize,
) !void {
    const abs_incx = @abs(incx);

    const x_expected = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x_expected);

    const x_actual = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x_actual);

    var alpha: N = tzsl.randomNumber(N, rand);

    for (x_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);

        x_actual[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dscal(
            zsl.numeric.cast(isize, n),
            alpha,
            x_expected.ptr,
            zsl.numeric.cast(isize, abs_incx),
        )
    else
        zsl.linalg.cblas.zscal(
            zsl.numeric.cast(isize, n),
            &alpha,
            x_expected.ptr,
            zsl.numeric.cast(isize, abs_incx),
        );

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
        std.debug.print("\n\tSCAL Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {}\n", .{ @typeName(N), n, incx });
        return e;
    };

    if (n == 0) return;

    const rel_tol = 1e-14;
    const abs_tol = 1e-10;
    for (x_expected, x_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tSCAL Test Failed\n", .{});
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
                std.debug.print("\n\tSCAL Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | index: {}\n", .{ n, incx, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:             {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
