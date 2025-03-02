const std = @import("std");
const core = @import("../../core/core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

pub inline fn syr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("blas.syr2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        A[@intCast(iaij)] += t0 * y[@intCast(iy)] + t1 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    A[@intCast(iaij)] += t0 * y[@intCast(jy)] + t1 * x[@intCast(jx)];

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            } else {
                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)] += t0 * y[@intCast(jy)] + t1 * x[@intCast(jx)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    var iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        A[@intCast(iaij)] += t0 * y[@intCast(iy)] + t1 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .Complex => @compileError("blas.syr2 does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.syr2 only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "syr2" {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, 2 * n);
    defer a.free(y1);

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
    @memcpy(y1.ptr, &[_]f64{
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

    blas.syr2(f64, .RowMajor, .Upper, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, n);

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

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
    defer a.free(y2);

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
    @memcpy(y2.ptr, &[_]f64{
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

    blas.syr2(f64, .RowMajor, .Upper, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(9, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(6, A[5]);
    try std.testing.expectEqual(39, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(90, A[9]);
    try std.testing.expectEqual(11, A[10]);
    try std.testing.expectEqual(12, A[11]);
    try std.testing.expectEqual(85, A[12]);
    try std.testing.expectEqual(110, A[13]);
    try std.testing.expectEqual(135, A[14]);
    try std.testing.expectEqual(16, A[15]);
    try std.testing.expectEqual(17, A[16]);
    try std.testing.expectEqual(18, A[17]);
    try std.testing.expectEqual(147, A[18]);
    try std.testing.expectEqual(180, A[19]);
    try std.testing.expectEqual(21, A[20]);
    try std.testing.expectEqual(22, A[21]);
    try std.testing.expectEqual(23, A[22]);
    try std.testing.expectEqual(24, A[23]);
    try std.testing.expectEqual(225, A[24]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

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
    @memcpy(y3.ptr, &[_]f64{
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

    blas.syr2(f64, .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(13, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(14, A[5]);
    try std.testing.expectEqual(55, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(90, A[9]);
    try std.testing.expectEqual(23, A[10]);
    try std.testing.expectEqual(36, A[11]);
    try std.testing.expectEqual(121, A[12]);
    try std.testing.expectEqual(110, A[13]);
    try std.testing.expectEqual(135, A[14]);
    try std.testing.expectEqual(32, A[15]);
    try std.testing.expectEqual(49, A[16]);
    try std.testing.expectEqual(66, A[17]);
    try std.testing.expectEqual(211, A[18]);
    try std.testing.expectEqual(180, A[19]);
    try std.testing.expectEqual(41, A[20]);
    try std.testing.expectEqual(62, A[21]);
    try std.testing.expectEqual(83, A[22]);
    try std.testing.expectEqual(104, A[23]);
    try std.testing.expectEqual(325, A[24]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * n);
    defer a.free(y4);

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
    @memcpy(y4.ptr, &[_]f64{
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

    blas.syr2(f64, .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, n);

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

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);
    const y5 = try a.alloc(f64, 2 * n);
    defer a.free(y5);

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
    @memcpy(y5.ptr, &[_]f64{
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

    blas.syr2(f64, .RowMajor, .Lower, n, alpha, x5.ptr, 2, y5.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(21, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(30, A[5]);
    try std.testing.expectEqual(87, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(90, A[9]);
    try std.testing.expectEqual(47, A[10]);
    try std.testing.expectEqual(84, A[11]);
    try std.testing.expectEqual(193, A[12]);
    try std.testing.expectEqual(110, A[13]);
    try std.testing.expectEqual(135, A[14]);
    try std.testing.expectEqual(64, A[15]);
    try std.testing.expectEqual(113, A[16]);
    try std.testing.expectEqual(162, A[17]);
    try std.testing.expectEqual(339, A[18]);
    try std.testing.expectEqual(180, A[19]);
    try std.testing.expectEqual(81, A[20]);
    try std.testing.expectEqual(142, A[21]);
    try std.testing.expectEqual(203, A[22]);
    try std.testing.expectEqual(264, A[23]);
    try std.testing.expectEqual(525, A[24]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);
    const y6 = try a.alloc(f64, 2 * n);
    defer a.free(y6);

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
    @memcpy(y6.ptr, &[_]f64{
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

    blas.syr2(f64, .RowMajor, .Lower, n, alpha, x6.ptr, -2, y6.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(25, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(38, A[5]);
    try std.testing.expectEqual(103, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(90, A[9]);
    try std.testing.expectEqual(59, A[10]);
    try std.testing.expectEqual(108, A[11]);
    try std.testing.expectEqual(229, A[12]);
    try std.testing.expectEqual(110, A[13]);
    try std.testing.expectEqual(135, A[14]);
    try std.testing.expectEqual(80, A[15]);
    try std.testing.expectEqual(145, A[16]);
    try std.testing.expectEqual(210, A[17]);
    try std.testing.expectEqual(403, A[18]);
    try std.testing.expectEqual(180, A[19]);
    try std.testing.expectEqual(101, A[20]);
    try std.testing.expectEqual(182, A[21]);
    try std.testing.expectEqual(263, A[22]);
    try std.testing.expectEqual(344, A[23]);
    try std.testing.expectEqual(625, A[24]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);
    const y7 = try a.alloc(f64, 2 * n);
    defer a.free(y7);

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
    @memcpy(y7.ptr, &[_]f64{
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

    blas.syr2(f64, .ColumnMajor, .Lower, n, alpha, x7.ptr, 2, y7.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(29, A[0]);
    try std.testing.expectEqual(26, A[1]);
    try std.testing.expectEqual(39, A[2]);
    try std.testing.expectEqual(52, A[3]);
    try std.testing.expectEqual(65, A[4]);
    try std.testing.expectEqual(38, A[5]);
    try std.testing.expectEqual(119, A[6]);
    try std.testing.expectEqual(80, A[7]);
    try std.testing.expectEqual(105, A[8]);
    try std.testing.expectEqual(130, A[9]);
    try std.testing.expectEqual(59, A[10]);
    try std.testing.expectEqual(108, A[11]);
    try std.testing.expectEqual(265, A[12]);
    try std.testing.expectEqual(158, A[13]);
    try std.testing.expectEqual(195, A[14]);
    try std.testing.expectEqual(80, A[15]);
    try std.testing.expectEqual(145, A[16]);
    try std.testing.expectEqual(210, A[17]);
    try std.testing.expectEqual(467, A[18]);
    try std.testing.expectEqual(260, A[19]);
    try std.testing.expectEqual(101, A[20]);
    try std.testing.expectEqual(182, A[21]);
    try std.testing.expectEqual(263, A[22]);
    try std.testing.expectEqual(344, A[23]);
    try std.testing.expectEqual(725, A[24]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);
    const y8 = try a.alloc(f64, 2 * n);
    defer a.free(y8);

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
    @memcpy(y8.ptr, &[_]f64{
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

    blas.syr2(f64, .ColumnMajor, .Lower, n, alpha, x8.ptr, -2, y8.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(33, A[0]);
    try std.testing.expectEqual(34, A[1]);
    try std.testing.expectEqual(51, A[2]);
    try std.testing.expectEqual(68, A[3]);
    try std.testing.expectEqual(85, A[4]);
    try std.testing.expectEqual(38, A[5]);
    try std.testing.expectEqual(135, A[6]);
    try std.testing.expectEqual(104, A[7]);
    try std.testing.expectEqual(137, A[8]);
    try std.testing.expectEqual(170, A[9]);
    try std.testing.expectEqual(59, A[10]);
    try std.testing.expectEqual(108, A[11]);
    try std.testing.expectEqual(301, A[12]);
    try std.testing.expectEqual(206, A[13]);
    try std.testing.expectEqual(255, A[14]);
    try std.testing.expectEqual(80, A[15]);
    try std.testing.expectEqual(145, A[16]);
    try std.testing.expectEqual(210, A[17]);
    try std.testing.expectEqual(531, A[18]);
    try std.testing.expectEqual(340, A[19]);
    try std.testing.expectEqual(101, A[20]);
    try std.testing.expectEqual(182, A[21]);
    try std.testing.expectEqual(263, A[22]);
    try std.testing.expectEqual(344, A[23]);
    try std.testing.expectEqual(825, A[24]);
}
