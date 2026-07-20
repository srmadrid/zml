const std = @import("std");

const zsl = @import("zsl");
const copy = zsl.linalg.blas.copy;
const copyParallel = zsl.linalg.blas.copyParallel;

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

test copy {
    @setEvalBranchQuota(10_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toSeconds()));
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeCopyTest(
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

fn executeCopyTest(
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

    const y_expected = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_expected);

    const y_actual = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_actual);
    const y_actual_parallel = try allocator.alloc(N, if (n == 0) 0 else 1 + (n - 1) * abs_incy);
    defer allocator.free(y_actual_parallel);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);

        y_actual[i] = elem.*;
        y_actual_parallel[i] = elem.*;
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
    ) catch |e| {
        std.debug.print("\n\tCOPY Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    copyParallel(
        n,
        x.ptr,
        incx,
        y_actual_parallel.ptr,
        incy,
        pool,
    ) catch |e| {
        std.debug.print("\n\tCOPY Test Failed\n", .{});
        std.debug.print("Type: {s} | n: {} | incx: {} | incy: {}\n", .{ @typeName(N), n, incx, incy });
        return e;
    };

    if (n == 0)
        return;

    for (y_expected, y_actual, y_actual_parallel, 0..) |expected, actual, actual_parallel, i| {
        if (N == f64) {
            if (expected != actual or expected != actual_parallel) {
                std.debug.print("\n\tCOPY Test Failed\n", .{});
                std.debug.print("Type: f64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (zsl):       {d}\n", .{actual});
                std.debug.print("Actual (zsl, pool): {d}\n", .{actual});

                return error.TestFailed;
            }
        } else {
            if (expected.re != actual.re or expected.im != actual.im or
                expected.re != actual_parallel.re or expected.im != actual_parallel.im)
            {
                std.debug.print("\n\tCOPY Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | n: {} | incx: {} | incy: {} | index: {}\n", .{ n, incx, incy, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Actual (zsl, pool): {d} + {d}i\n", .{ actual_parallel.re, actual_parallel.im });

                return error.TestFailed;
            }
        }
    }
}
