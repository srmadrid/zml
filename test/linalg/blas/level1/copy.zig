const std = @import("std");

const zsl = @import("zsl");
const copy = zsl.linalg.blas.copy;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ @as(usize, 0), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 0), @as(isize, -1), @as(isize, 2), null, null },

    // Single element
    .{ @as(usize, 1), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 1), @as(isize, 2), @as(isize, -1), null, null },
    .{ @as(usize, 1), @as(isize, -1), @as(isize, 2), null, null },

    // unroll aligned
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 64), @as(isize, 1), @as(isize, 1), null, null },

    // unroll unaligned / remainder loops
    .{ @as(usize, 37), @as(isize, 1), @as(isize, 1), null, null },
    .{ @as(usize, 69), @as(isize, 1), @as(isize, 1), null, null },

    // Strided
    .{ @as(usize, 32), @as(isize, 2), @as(isize, 1), null, null },
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 2), null, null },
    .{ @as(usize, 65), @as(isize, 3), @as(isize, 3), null, null },
    .{ @as(usize, 32), @as(isize, -1), @as(isize, 1), null, null },
    .{ @as(usize, 65), @as(isize, 1), @as(isize, -2), null, null },
    .{ @as(usize, 65), @as(isize, -2), @as(isize, -1), null, null },

    // Lowering parallel_threshold to test threading logic on small data chunks
    .{ @as(usize, 32), @as(isize, 1), @as(isize, 1), null, @as(usize, 10) },
    .{ @as(usize, 33), @as(isize, 1), @as(isize, 1), null, @as(usize, 8) },

    // Forced multithreaded
    .{ @as(usize, 33), @as(isize, 2), @as(isize, -1), @as(usize, 2), null },
    .{ @as(usize, 33), @as(isize, -1), @as(isize, 2), @as(usize, 4), null },
    .{ @as(usize, 33), @as(isize, -3), @as(isize, -3), @as(usize, 2), null },

    // Default threading behavior
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

test copy {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeCopyTest(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                combo[3],
                combo[4],
            );
        }
    }
}

fn executeCopyTest(
    comptime N: type,
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

    const x = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incx);
    defer allocator.free(x);

    const y_expected = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_expected);

    const y_actual = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_actual);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);

        y_actual[i] = elem.*;
    }

    if (N == f64)
        zsl.linalg.cblas.dcopy(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
        )
    else
        zsl.linalg.cblas.zcopy(
            zsl.numeric.cast(isize, n),
            x.ptr,
            zsl.numeric.cast(isize, incx),
            y_expected.ptr,
            zsl.numeric.cast(isize, incy),
        );

    copy(
        n,
        x.ptr,
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
        std.debug.print("\n\tCOPY Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    if (n == 0) return;

    for (y_expected, y_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            if (expected != actual) {
                std.debug.print("\n\tCOPY Test Failed\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d}\n", .{expected});
                std.debug.print("Actual (zsl):     {d}\n", .{actual});

                return error.TestFailed;
            }
        } else {
            if (expected.re != actual.re or expected.im != actual.im) {
                std.debug.print("\n\tCOPY Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS): {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):     {d} + {d}i\n", .{ actual.re, actual.im });

                return error.TestFailed;
            }
        }
    }
}
