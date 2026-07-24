const std = @import("std");

const zsl = @import("zsl");
const syr2 = zsl.linalg.blas.syr2;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays
    .{ .col_major, .upper, @as(usize, 0), @as(isize, 1), @as(isize, 1), @as(usize, 1) },
    .{ .row_major, .lower, @as(usize, 0), @as(isize, -1), @as(isize, 2), @as(usize, 10) },

    // Single element
    .{ .col_major, .upper, @as(usize, 1), @as(isize, 1), @as(isize, 1), @as(usize, 1) },
    .{ .row_major, .lower, @as(usize, 1), @as(isize, -1), @as(isize, -1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .upper, @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 32) },
    .{ .row_major, .lower, @as(usize, 64), @as(isize, 1), @as(isize, 1), @as(usize, 64) },

    // Requires lda >= n
    .{ .col_major, .upper, @as(usize, 32), @as(isize, 1), @as(isize, 1), @as(usize, 40) },
    .{ .row_major, .lower, @as(usize, 37), @as(isize, 1), @as(isize, 1), @as(usize, 45) },

    // Strided
    .{ .col_major, .lower, @as(usize, 32), @as(isize, 2), @as(isize, 1), @as(usize, 32) },
    .{ .row_major, .upper, @as(usize, 32), @as(isize, 3), @as(isize, 2), @as(usize, 32) },
    .{ .col_major, .upper, @as(usize, 37), @as(isize, -1), @as(isize, 2), @as(usize, 40) },
    .{ .row_major, .lower, @as(usize, 37), @as(isize, -2), @as(isize, -1), @as(usize, 40) },
    .{ .col_major, .lower, @as(usize, 37), @as(isize, -2), @as(isize, -2), @as(usize, 37) },
    .{ .col_major, .upper, @as(usize, 64), @as(isize, 5), @as(isize, 1), @as(usize, 64) },
    .{ .row_major, .lower, @as(usize, 65), @as(isize, 1), @as(isize, -1), @as(usize, 65) },
    .{ .col_major, .lower, @as(usize, 65), @as(isize, 2), @as(isize, 2), @as(usize, 65) },
    .{ .row_major, .upper, @as(usize, 65), @as(isize, -1), @as(isize, 3), @as(usize, 70) },
    .{ .col_major, .upper, @as(usize, 65), @as(isize, -3), @as(isize, -2), @as(usize, 65) },

    // Large arrays
    .{ .col_major, .upper, @as(usize, 2000), @as(isize, 1), @as(isize, 1), @as(usize, 2000) },
    .{ .row_major, .lower, @as(usize, 1500), @as(isize, 2), @as(isize, -1), @as(usize, 1500) },
    .{ .row_major, .lower, @as(usize, 1500), @as(isize, 1), @as(isize, 1), @as(usize, 1500) },
    .{ .col_major, .upper, @as(usize, 1000), @as(isize, 3), @as(isize, 2), @as(usize, 1005) },
    .{ .col_major, .lower, @as(usize, 1500), @as(isize, 1), @as(isize, 1), @as(usize, 1500) },
    .{ .row_major, .upper, @as(usize, 1000), @as(isize, -1), @as(isize, 1), @as(usize, 1000) },

    // Prime-like numbers to ensure threads get unequal chunks
    .{ .col_major, .upper, @as(usize, 1507), @as(isize, 1), @as(isize, 1), @as(usize, 1510) },
    .{ .row_major, .lower, @as(usize, 1003), @as(isize, 2), @as(isize, -2), @as(usize, 1010) },
};

test syr2 {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    var pool: zsl.thread.Pool = undefined;
    try pool.init(allocator, .{});
    defer pool.deinit(allocator);

    inline for (combinations) |combo| {
        try executeSyr2Test(
            f64,
            allocator,
            rand,
            combo[0],
            combo[1],
            combo[2],
            combo[3],
            combo[4],
            combo[5],
            &pool,
        );
    }
}

fn executeSyr2Test(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    uplo: zsl.matrix.Uplo,
    n: usize,
    incx: isize,
    incy: isize,
    lda: usize,
    pool: *zsl.thread.Pool,
) !void {
    _ = pool;

    const abs_incx = @abs(incx);
    const abs_incy = @abs(incy);

    const len_mem_x = if (n == 0) 0 else 1 + (n - 1) * abs_incx;
    const len_mem_y = if (n == 0) 0 else 1 + (n - 1) * abs_incy;
    const len_mem_a = lda * n;

    const x = try allocator.alloc(N, len_mem_x);
    defer allocator.free(x);

    const y = try allocator.alloc(N, len_mem_y);
    defer allocator.free(y);

    const a_expected = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_expected);

    const a_actual = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a_actual);

    const alpha: N = tzsl.randomNumber(N, rand);

    for (x) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (y) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (a_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        a_actual[i] = elem.*;
    }

    zsl.linalg.cblas.dsyr2(
        layout.toInt(c_int),
        uplo.toInt(c_int),
        zsl.numeric.cast(isize, n),
        alpha,
        x.ptr,
        zsl.numeric.cast(isize, incx),
        y.ptr,
        zsl.numeric.cast(isize, incy),
        a_expected.ptr,
        zsl.numeric.cast(isize, lda),
    );

    syr2(
        layout,
        uplo,
        n,
        alpha,
        x.ptr,
        incx,
        y.ptr,
        incy,
        a_actual.ptr,
        lda,
    ) catch |e| {
        std.debug.print("\n\tSYR2 Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | n: {} | incx: {} | incy: {} | lda: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), n, incx, incy, lda });
        return e;
    };

    if (n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, n);
    const abs_tol = 1e-10;
    for (a_expected, a_actual, 0..) |expected, actual, i| {
        const diff = zsl.float.abs(expected - actual);
        if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
            std.debug.print("\n\tSYR2 Test Failed\n", .{});
            std.debug.print("Type: f64 | layout: {s} | uplo: {s} | n: {} | incx: {} | lda: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), n, incx, lda, i });
            std.debug.print("Expected (CBLAS): {d}\n", .{expected});
            std.debug.print("Actual (zsl):     {d}\n", .{actual});
            std.debug.print("Diff:             {d}\n", .{diff});

            return error.TestFailed;
        }
    }
}
