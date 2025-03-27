const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

pub inline fn syr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.syr does not support bool."),
        .int, .float => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < j) {
                        A[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    A[@intCast(iaij)] += t0 * x[@intCast(jx)];

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)] += t0 * x[@intCast(jx)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        A[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                }
            }
        },
        .cfloat => @compileError("blas.syr does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.syr only supports simple types."),
        .unsupported => unreachable,
    }
}

test syr {
    const a = std.testing.allocator;

    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);

    @memcpy(A.ptr, &[_]f64{
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    });
    @memcpy(x1.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    blas.syr(f64, .RowMajor, .Upper, n, alpha, x1.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(3, A[0]);
    try std.testing.expectEqual(6, A[1]);
    try std.testing.expectEqual(9, A[2]);
    try std.testing.expectEqual(12, A[3]);
    try std.testing.expectEqual(15, A[4]);
    try std.testing.expectEqual(6, A[5]);
    try std.testing.expectEqual(15, A[6]);
    try std.testing.expectEqual(20, A[7]);
    try std.testing.expectEqual(25, A[8]);
    try std.testing.expectEqual(30, A[9]);
    try std.testing.expectEqual(11, A[10]);
    try std.testing.expectEqual(12, A[11]);
    try std.testing.expectEqual(31, A[12]);
    try std.testing.expectEqual(38, A[13]);
    try std.testing.expectEqual(45, A[14]);
    try std.testing.expectEqual(16, A[15]);
    try std.testing.expectEqual(17, A[16]);
    try std.testing.expectEqual(18, A[17]);
    try std.testing.expectEqual(51, A[18]);
    try std.testing.expectEqual(60, A[19]);
    try std.testing.expectEqual(21, A[20]);
    try std.testing.expectEqual(22, A[21]);
    try std.testing.expectEqual(23, A[22]);
    try std.testing.expectEqual(24, A[23]);
    try std.testing.expectEqual(75, A[24]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.syr(f64, .RowMajor, .Upper, n, alpha, x2.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(5, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(6, A[5]);
    try std.testing.expectEqual(23, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(11, A[10]);
    try std.testing.expectEqual(12, A[11]);
    try std.testing.expectEqual(49, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(16, A[15]);
    try std.testing.expectEqual(17, A[16]);
    try std.testing.expectEqual(18, A[17]);
    try std.testing.expectEqual(83, A[18]);
    try std.testing.expectEqual(100, A[19]);
    try std.testing.expectEqual(21, A[20]);
    try std.testing.expectEqual(22, A[21]);
    try std.testing.expectEqual(23, A[22]);
    try std.testing.expectEqual(24, A[23]);
    try std.testing.expectEqual(125, A[24]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    blas.syr(f64, .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(7, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(10, A[5]);
    try std.testing.expectEqual(31, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(17, A[10]);
    try std.testing.expectEqual(24, A[11]);
    try std.testing.expectEqual(67, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(24, A[15]);
    try std.testing.expectEqual(33, A[16]);
    try std.testing.expectEqual(42, A[17]);
    try std.testing.expectEqual(115, A[18]);
    try std.testing.expectEqual(100, A[19]);
    try std.testing.expectEqual(31, A[20]);
    try std.testing.expectEqual(42, A[21]);
    try std.testing.expectEqual(53, A[22]);
    try std.testing.expectEqual(64, A[23]);
    try std.testing.expectEqual(175, A[24]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);

    @memcpy(x4.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.syr(f64, .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(9, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(14, A[5]);
    try std.testing.expectEqual(39, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(23, A[10]);
    try std.testing.expectEqual(36, A[11]);
    try std.testing.expectEqual(85, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(32, A[15]);
    try std.testing.expectEqual(49, A[16]);
    try std.testing.expectEqual(66, A[17]);
    try std.testing.expectEqual(147, A[18]);
    try std.testing.expectEqual(100, A[19]);
    try std.testing.expectEqual(41, A[20]);
    try std.testing.expectEqual(62, A[21]);
    try std.testing.expectEqual(83, A[22]);
    try std.testing.expectEqual(104, A[23]);
    try std.testing.expectEqual(225, A[24]);

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    blas.syr(f64, .RowMajor, .Lower, n, alpha, x5.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(11, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(18, A[5]);
    try std.testing.expectEqual(47, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(29, A[10]);
    try std.testing.expectEqual(48, A[11]);
    try std.testing.expectEqual(103, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(40, A[15]);
    try std.testing.expectEqual(65, A[16]);
    try std.testing.expectEqual(90, A[17]);
    try std.testing.expectEqual(179, A[18]);
    try std.testing.expectEqual(100, A[19]);
    try std.testing.expectEqual(51, A[20]);
    try std.testing.expectEqual(82, A[21]);
    try std.testing.expectEqual(113, A[22]);
    try std.testing.expectEqual(144, A[23]);
    try std.testing.expectEqual(275, A[24]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.syr(f64, .RowMajor, .Lower, n, alpha, x6.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(13, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(22, A[5]);
    try std.testing.expectEqual(55, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(35, A[10]);
    try std.testing.expectEqual(60, A[11]);
    try std.testing.expectEqual(121, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(48, A[15]);
    try std.testing.expectEqual(81, A[16]);
    try std.testing.expectEqual(114, A[17]);
    try std.testing.expectEqual(211, A[18]);
    try std.testing.expectEqual(100, A[19]);
    try std.testing.expectEqual(61, A[20]);
    try std.testing.expectEqual(102, A[21]);
    try std.testing.expectEqual(143, A[22]);
    try std.testing.expectEqual(184, A[23]);
    try std.testing.expectEqual(325, A[24]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    blas.syr(f64, .ColumnMajor, .Lower, n, alpha, x7.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(15, A[0]);
    try std.testing.expectEqual(14, A[1]);
    try std.testing.expectEqual(21, A[2]);
    try std.testing.expectEqual(28, A[3]);
    try std.testing.expectEqual(35, A[4]);
    try std.testing.expectEqual(22, A[5]);
    try std.testing.expectEqual(63, A[6]);
    try std.testing.expectEqual(44, A[7]);
    try std.testing.expectEqual(57, A[8]);
    try std.testing.expectEqual(70, A[9]);
    try std.testing.expectEqual(35, A[10]);
    try std.testing.expectEqual(60, A[11]);
    try std.testing.expectEqual(139, A[12]);
    try std.testing.expectEqual(86, A[13]);
    try std.testing.expectEqual(105, A[14]);
    try std.testing.expectEqual(48, A[15]);
    try std.testing.expectEqual(81, A[16]);
    try std.testing.expectEqual(114, A[17]);
    try std.testing.expectEqual(243, A[18]);
    try std.testing.expectEqual(140, A[19]);
    try std.testing.expectEqual(61, A[20]);
    try std.testing.expectEqual(102, A[21]);
    try std.testing.expectEqual(143, A[22]);
    try std.testing.expectEqual(184, A[23]);
    try std.testing.expectEqual(375, A[24]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.syr(f64, .ColumnMajor, .Lower, n, alpha, x8.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(17, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(22, A[5]);
    try std.testing.expectEqual(71, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(90, A[9]);
    try std.testing.expectEqual(35, A[10]);
    try std.testing.expectEqual(60, A[11]);
    try std.testing.expectEqual(157, A[12]);
    try std.testing.expectEqual(110, A[13]);
    try std.testing.expectEqual(135, A[14]);
    try std.testing.expectEqual(48, A[15]);
    try std.testing.expectEqual(81, A[16]);
    try std.testing.expectEqual(114, A[17]);
    try std.testing.expectEqual(275, A[18]);
    try std.testing.expectEqual(180, A[19]);
    try std.testing.expectEqual(61, A[20]);
    try std.testing.expectEqual(102, A[21]);
    try std.testing.expectEqual(143, A[22]);
    try std.testing.expectEqual(184, A[23]);
    try std.testing.expectEqual(425, A[24]);
}
