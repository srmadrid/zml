const std = @import("std");

const zsl = @import("zsl");
const syrk = zsl.linalg.blas.syrk;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays (all zero)
    .{ .col_major, .upper, .trans, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .lower, .no_trans, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1) },

    // k = 0 (Testing diagonal/beta scaling fallback)
    .{ .col_major, .lower, .trans, @as(usize, 10), @as(usize, 0), @as(usize, 10), @as(usize, 10) },
    .{ .row_major, .upper, .trans, @as(usize, 12), @as(usize, 0), @as(usize, 12), @as(usize, 12) },

    // n = 0
    .{ .col_major, .upper, .no_trans, @as(usize, 0), @as(usize, 5), @as(usize, 1), @as(usize, 1) },

    // Single element
    .{ .col_major, .upper, .no_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .lower, .no_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .lower, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32) },
    .{ .row_major, .upper, .trans, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64) },

    // Rectangular / Transposed variants (A is n-by-k if no_trans, k-by-n if trans)
    .{ .col_major, .upper, .no_trans, @as(usize, 15), @as(usize, 25), @as(usize, 15), @as(usize, 15) },
    .{ .col_major, .lower, .trans, @as(usize, 15), @as(usize, 25), @as(usize, 25), @as(usize, 15) },
    .{ .row_major, .upper, .trans, @as(usize, 17), @as(usize, 13), @as(usize, 17), @as(usize, 17) },
    .{ .row_major, .lower, .no_trans, @as(usize, 17), @as(usize, 13), @as(usize, 17), @as(usize, 17) },

    // Padded leading dimensions (lda, ldc > minimum required)
    .{ .col_major, .upper, .no_trans, @as(usize, 10), @as(usize, 12), @as(usize, 15), @as(usize, 14) },
    .{ .row_major, .lower, .trans, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 15) },

    // Large arrays
    .{ .col_major, .upper, .no_trans, @as(usize, 100), @as(usize, 120), @as(usize, 100), @as(usize, 100) },
    .{ .row_major, .lower, .trans, @as(usize, 100), @as(usize, 120), @as(usize, 120), @as(usize, 100) },
    .{ .row_major, .lower, .trans, @as(usize, 1000), @as(usize, 1200), @as(usize, 1200), @as(usize, 1000) },

    // Prime-like numbers
    .{ .col_major, .lower, .trans, @as(usize, 67), @as(usize, 73), @as(usize, 75), @as(usize, 80) },
    .{ .row_major, .upper, .no_trans, @as(usize, 71), @as(usize, 67), @as(usize, 75), @as(usize, 75) },
};

test syrk {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeSyrkTest(
                N,
                allocator,
                rand,
                combo[0],
                combo[1],
                combo[2],
                combo[3],
                combo[4],
                combo[5],
                combo[6],
            );
        }
    }
}

fn executeSyrkTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    uplo: zsl.matrix.Uplo,
    transa: zsl.linalg.blas.Transpose,
    n: usize,
    k: usize,
    lda: usize,
    ldc: usize,
) !void {
    const notrans = transa == .no_trans or transa == .no_trans;
    const nrowa: usize = if (notrans) n else k;
    const ncola: usize = if (notrans) k else n;

    const len_mem_a = lda * if (layout == .col_major) ncola else nrowa;

    const len_mem_c = ldc * n;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const c_expected = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_expected);

    const c_actual = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_actual);

    const alpha: N = tzsl.randomNumber(N, rand);
    const beta: N = tzsl.randomNumber(N, rand);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (c_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        c_actual[i] = elem.*;
    }

    if (N == f64) {
        zsl.linalg.cblas.dsyrk(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );
    } else {
        zsl.linalg.cblas.zsyrk(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transa.toInt(c_int),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            &beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );
    }

    syrk(
        layout,
        uplo,
        transa,
        n,
        k,
        alpha,
        a.ptr,
        lda,
        beta,
        c_actual.ptr,
        ldc,
    ) catch |e| {
        std.debug.print("\n\tSYRK Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | transa: {s} | n: {} | k: {} | lda: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), @tagName(transa), n, k, lda, ldc });
        return e;
    };

    if (n == 0)
        return;

    const rel_tol = 1e-14 * zsl.numeric.cast(f64, zsl.int.max(1, k));
    const abs_tol = 1e-10;
    for (c_expected, c_actual, 0..) |expected, actual, i| {
        if (N == f64) {
            const diff = zsl.float.abs(expected - actual);
            if (diff > abs_tol and diff > zsl.float.abs(expected) * rel_tol) {
                std.debug.print("\n\tSYRK Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | uplo: {s} | transa: {s} | n: {} | k: {} | lda: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transa), n, k, lda, ldc, i });
                std.debug.print("Expected (CBLAS):   {d}\n", .{expected});
                std.debug.print("Actual (zsl):       {d}\n", .{actual});
                std.debug.print("Diff:               {d}\n", .{diff});

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
                std.debug.print("\n\tSYRK Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | uplo: {s} | transa: {s} | n: {} | k: {} | lda: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transa), n, k, lda, ldc, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
