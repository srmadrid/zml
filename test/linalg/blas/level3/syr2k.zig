const std = @import("std");

const zsl = @import("zsl");
const syr2k = zsl.linalg.blas.syr2k;

const tzsl = @import("../../../zsl.zig");

const combinations = .{
    // Empty arrays (all zero)
    .{ .col_major, .upper, .no_trans, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .lower, .trans, @as(usize, 0), @as(usize, 0), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // k = 0 (Testing diagonal/beta scaling fallback)
    .{ .col_major, .lower, .no_trans, @as(usize, 10), @as(usize, 0), @as(usize, 10), @as(usize, 10), @as(usize, 10) },
    .{ .row_major, .upper, .trans, @as(usize, 12), @as(usize, 0), @as(usize, 12), @as(usize, 12), @as(usize, 12) },

    // n = 0
    .{ .col_major, .upper, .no_trans, @as(usize, 0), @as(usize, 5), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Single element
    .{ .col_major, .upper, .trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },
    .{ .row_major, .lower, .no_trans, @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1), @as(usize, 1) },

    // Square aligned
    .{ .col_major, .lower, .no_trans, @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32), @as(usize, 32) },
    .{ .row_major, .upper, .trans, @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64), @as(usize, 64) },

    // Rectangular / Transposed variants (A and B are n-by-k if no_trans, k-by-n if trans)
    .{ .col_major, .upper, .no_trans, @as(usize, 15), @as(usize, 25), @as(usize, 15), @as(usize, 15), @as(usize, 15) },
    .{ .col_major, .lower, .trans, @as(usize, 15), @as(usize, 25), @as(usize, 25), @as(usize, 25), @as(usize, 15) },
    .{ .row_major, .upper, .trans, @as(usize, 17), @as(usize, 13), @as(usize, 17), @as(usize, 17), @as(usize, 17) },
    .{ .row_major, .lower, .no_trans, @as(usize, 17), @as(usize, 13), @as(usize, 17), @as(usize, 17), @as(usize, 17) },

    // Padded leading dimensions (lda, ldb, ldc > minimum required; lda != ldb to catch swap bugs)
    .{ .col_major, .upper, .no_trans, @as(usize, 10), @as(usize, 12), @as(usize, 15), @as(usize, 16), @as(usize, 14) },
    .{ .row_major, .lower, .trans, @as(usize, 10), @as(usize, 12), @as(usize, 14), @as(usize, 13), @as(usize, 15) },

    // Large arrays
    .{ .col_major, .upper, .no_trans, @as(usize, 100), @as(usize, 120), @as(usize, 100), @as(usize, 100), @as(usize, 100) },
    .{ .row_major, .lower, .trans, @as(usize, 100), @as(usize, 120), @as(usize, 120), @as(usize, 120), @as(usize, 100) },
    .{ .row_major, .lower, .no_trans, @as(usize, 1000), @as(usize, 1200), @as(usize, 1200), @as(usize, 1200), @as(usize, 1000) },

    // Prime-like numbers
    .{ .col_major, .lower, .trans, @as(usize, 67), @as(usize, 73), @as(usize, 87), @as(usize, 89), @as(usize, 77) },
    .{ .row_major, .upper, .no_trans, @as(usize, 71), @as(usize, 67), @as(usize, 75), @as(usize, 77), @as(usize, 75) },
};

test syr2k {
    @setEvalBranchQuota(20_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    inline for (combinations) |combo| {
        inline for (.{ f64, zsl.cf64 }) |N| {
            try executeSyr2kTest(
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
                combo[7],
            );
        }
    }
}

fn executeSyr2kTest(
    comptime N: type,
    allocator: std.mem.Allocator,
    rand: std.Random,
    layout: zsl.matrix.Layout,
    uplo: zsl.matrix.Uplo,
    transab: zsl.linalg.blas.Transpose,
    n: usize,
    k: usize,
    lda: usize,
    ldb: usize,
    ldc: usize,
) !void {
    const notrans = transab == .no_trans or transab == .conj_no_trans;
    const nrowab: usize = if (notrans) n else k;
    const ncolab: usize = if (notrans) k else n;

    const len_mem_a = lda * if (layout == .col_major) ncolab else nrowab;
    const len_mem_b = ldb * if (layout == .col_major) ncolab else nrowab;
    const len_mem_c = ldc * n;

    const a = try allocator.alloc(N, len_mem_a);
    defer allocator.free(a);

    const b = try allocator.alloc(N, len_mem_b);
    defer allocator.free(b);

    const c_expected = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_expected);

    const c_actual = try allocator.alloc(N, len_mem_c);
    defer allocator.free(c_actual);

    const alpha: N = tzsl.randomNumber(N, rand);
    const beta: N = tzsl.randomNumber(N, rand);

    for (a) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (b) |*elem|
        elem.* = tzsl.randomNumber(N, rand);

    for (c_expected, 0..) |*elem, i| {
        elem.* = tzsl.randomNumber(N, rand);
        c_actual[i] = elem.*;
    }

    if (N == f64) {
        zsl.linalg.cblas.dsyr2k(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transab.toInt(c_int),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b.ptr,
            zsl.numeric.cast(isize, ldb),
            beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );
    } else {
        zsl.linalg.cblas.zsyr2k(
            layout.toInt(c_int),
            uplo.toInt(c_int),
            transab.toInt(c_int),
            zsl.numeric.cast(isize, n),
            zsl.numeric.cast(isize, k),
            &alpha,
            a.ptr,
            zsl.numeric.cast(isize, lda),
            b.ptr,
            zsl.numeric.cast(isize, ldb),
            &beta,
            c_expected.ptr,
            zsl.numeric.cast(isize, ldc),
        );
    }

    syr2k(
        layout,
        uplo,
        transab,
        n,
        k,
        alpha,
        a.ptr,
        lda,
        b.ptr,
        ldb,
        beta,
        c_actual.ptr,
        ldc,
    ) catch |e| {
        std.debug.print("\n\tSYR2K Test Failed\n", .{});
        std.debug.print("Type: {s} | layout: {s} | uplo: {s} | transab: {s} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {}\n", .{ @typeName(N), @tagName(layout), @tagName(uplo), @tagName(transab), n, k, lda, ldb, ldc });
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
                std.debug.print("\n\tSYR2K Test Failed\n", .{});
                std.debug.print("Type: f64 | layout: {s} | uplo: {s} | transab: {s} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transab), n, k, lda, ldb, ldc, i });
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
                std.debug.print("\n\tSYR2K Test Failed\n", .{});
                std.debug.print("Type: zsl.cf64 | layout: {s} | uplo: {s} | transab: {s} | n: {} | k: {} | lda: {} | ldb: {} | ldc: {} | index: {}\n", .{ @tagName(layout), @tagName(uplo), @tagName(transab), n, k, lda, ldb, ldc, i });
                std.debug.print("Expected (CBLAS):   {d} + {d}i\n", .{ expected.re, expected.im });
                std.debug.print("Actual (zsl):       {d} + {d}i\n", .{ actual.re, actual.im });
                std.debug.print("Diff:               {d} + {d}i\n", .{ diff_re, diff_im });

                return error.TestFailed;
            }
        }
    }
}
