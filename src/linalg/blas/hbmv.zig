const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Numeric = core.types.Numeric;

pub inline fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or k < 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: Numeric(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < k + 1) return;

    const LENX = N;
    const LENY = N;

    switch (numericType) {
        .bool => @compileError("blas.hbmv does not support bool."),
        .int, .float => @compileError("blas.hbmv does not support int or float."),
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0 and beta.re == 1 and beta.im == 0) return;

            if (alpha.re == 0 and alpha.im == 0) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta.re == 0 and beta.im == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)].re = 0;
                        y[@intCast(iy)].im = 0;

                        iy += incy;
                    }
                } else if (beta.re != 1 or beta.im != 0) {
                    while (iy != Sty) {
                        const tmp = y[@intCast(iy)].re * beta.re - y[@intCast(iy)].im * beta.im;
                        y[@intCast(iy)].im = y[@intCast(iy)].re * beta.im + y[@intCast(iy)].im * beta.re;
                        y[@intCast(iy)].re = tmp;

                        iy += incy;
                    }
                }

                return;
            }

            if (UPLO == .Upper) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta.re == 0 and beta.im == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)].re = 0;
                        y[@intCast(iy)].im = 0;

                        iy += incy;
                    }
                } else if (beta.re != 1 or beta.im != 0) {
                    while (iy != Sty) {
                        const tmp = y[@intCast(iy)].re * beta.re - y[@intCast(iy)].im * beta.im;
                        y[@intCast(iy)].im = y[@intCast(iy)].re * beta.im + y[@intCast(iy)].im * beta.re;
                        y[@intCast(iy)].re = tmp;

                        iy += incy;
                    }
                }

                var j: isize = 0;
                var jaj: isize = k;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    var t1 = T.init(0, 0);
                    const I0: isize = if (j - k < 0) 0 else j - k;
                    const I1: isize = j;

                    var i: isize = I0;
                    var iaij: isize = jaj - j + I0;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im * conj;
                        y[@intCast(iy)].im += t0.im * A[@intCast(iaij)].re + t0.re * A[@intCast(iaij)].im * conj;

                        t1.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * conj * x[@intCast(ix)].im;
                        t1.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * conj * x[@intCast(ix)].re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)].re += t0.re * A[@intCast(jaj)].re + alpha.re * t1.re - alpha.im * t1.im;
                    y[@intCast(jy)].im += t0.im * A[@intCast(jaj)].re + alpha.re * t1.im + alpha.im * t1.re;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            } else {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta.re == 0 and beta.im == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)].re = 0;
                        y[@intCast(iy)].im = 0;

                        iy += incy;
                    }
                } else if (beta.re != 1 or beta.im != 0) {
                    while (iy != Sty) {
                        const tmp = y[@intCast(iy)].re * beta.re - y[@intCast(iy)].im * beta.im;
                        y[@intCast(iy)].im = y[@intCast(iy)].re * beta.im + y[@intCast(iy)].im * beta.re;
                        y[@intCast(iy)].re = tmp;

                        iy += incy;
                    }
                }

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    var t1 = T.init(0, 0);
                    const I0: isize = j + 1;
                    const I1: isize = if (N - 1 > j + k) j + k + 1 else N;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im * conj;
                        y[@intCast(iy)].im += t0.im * A[@intCast(iaij)].re + t0.re * A[@intCast(iaij)].im * conj;

                        t1.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * conj * x[@intCast(ix)].im;
                        t1.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * conj * x[@intCast(ix)].re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)].re += t0.re * A[@intCast(jaj)].re + alpha.re * t1.re - alpha.im * t1.im;
                    y[@intCast(jy)].im += t0.im * A[@intCast(jaj)].re + alpha.re * t1.im + alpha.im * t1.re;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.hbmv only supports simple types."),
        .unsupported => unreachable,
    }
}

test "hbmv" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 1;
    const alpha = Complex(f64).init(1, 1);
    const beta = Complex(f64).init(3, 3);

    const A = try a.alloc(Complex(f64), (1 + 2 * k) * n);
    defer a.free(A);
    const x1 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(2, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(4, 4),
        Complex(f64).init(5, 5),
        Complex(f64).init(6, 6),
        Complex(f64).init(7, 7),
        Complex(f64).init(8, 8),
        Complex(f64).init(9, 9),
        Complex(f64).init(10, 10),
        Complex(f64).init(11, 11),
        Complex(f64).init(12, 12),
        Complex(f64).init(13, 13),
        Complex(f64).init(14, 14),
        Complex(f64).init(15, 15),
    });
    @memcpy(x1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });
    @memcpy(y1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .RowMajor, .Upper, n, k, alpha, A.ptr, k + 1, x1.ptr, 2, beta, y1.ptr, 2);

    try std.testing.expectEqual(-17, y1[0].re);
    try std.testing.expectEqual(21, y1[0].im);
    try std.testing.expectEqual(-47, y1[2].re);
    try std.testing.expectEqual(81, y1[2].im);
    try std.testing.expectEqual(-77, y1[4].re);
    try std.testing.expectEqual(189, y1[4].im);
    try std.testing.expectEqual(-107, y1[6].re);
    try std.testing.expectEqual(345, y1[6].im);
    try std.testing.expectEqual(103, y1[8].re);
    try std.testing.expectEqual(329, y1[8].im);

    const x2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y2.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .RowMajor, .Upper, n, k, alpha, A.ptr, k + 1, x2.ptr, -2, beta, y2.ptr, -2);

    try std.testing.expectEqual(-17, y2[8].re);
    try std.testing.expectEqual(21, y2[8].im);
    try std.testing.expectEqual(-47, y2[6].re);
    try std.testing.expectEqual(81, y2[6].im);
    try std.testing.expectEqual(-77, y2[4].re);
    try std.testing.expectEqual(189, y2[4].im);
    try std.testing.expectEqual(-107, y2[2].re);
    try std.testing.expectEqual(345, y2[2].im);
    try std.testing.expectEqual(103, y2[0].re);
    try std.testing.expectEqual(329, y2[0].im);

    const x3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });
    @memcpy(y3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .ColumnMajor, .Upper, n, k, alpha, A.ptr, k + 1, x3.ptr, 2, beta, y3.ptr, 2);

    try std.testing.expectEqual(-26, y3[0].re);
    try std.testing.expectEqual(30, y3[0].im);
    try std.testing.expectEqual(-58, y3[2].re);
    try std.testing.expectEqual(102, y3[2].im);
    try std.testing.expectEqual(-88, y3[4].re);
    try std.testing.expectEqual(222, y3[4].im);
    try std.testing.expectEqual(-118, y3[6].re);
    try std.testing.expectEqual(390, y3[6].im);
    try std.testing.expectEqual(116, y3[8].re);
    try std.testing.expectEqual(364, y3[8].im);

    const x4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y4.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .ColumnMajor, .Upper, n, k, alpha, A.ptr, k + 1, x4.ptr, -2, beta, y4.ptr, -2);

    try std.testing.expectEqual(-26, y4[8].re);
    try std.testing.expectEqual(30, y4[8].im);
    try std.testing.expectEqual(-58, y4[6].re);
    try std.testing.expectEqual(102, y4[6].im);
    try std.testing.expectEqual(-88, y4[4].re);
    try std.testing.expectEqual(222, y4[4].im);
    try std.testing.expectEqual(-118, y4[2].re);
    try std.testing.expectEqual(390, y4[2].im);
    try std.testing.expectEqual(116, y4[0].re);
    try std.testing.expectEqual(364, y4[0].im);

    const x5 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x5);
    const y5 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y5);

    @memcpy(x5.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });
    @memcpy(y5.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .RowMajor, .Lower, n, k, alpha, A.ptr, k + 1, x5.ptr, 2, beta, y5.ptr, 2);

    try std.testing.expectEqual(16, y5[0].re);
    try std.testing.expectEqual(36, y5[0].im);
    try std.testing.expectEqual(34, y5[2].re);
    try std.testing.expectEqual(106, y5[2].im);
    try std.testing.expectEqual(52, y5[4].re);
    try std.testing.expectEqual(226, y5[4].im);
    try std.testing.expectEqual(70, y5[6].re);
    try std.testing.expectEqual(394, y5[6].im);
    try std.testing.expectEqual(-154, y5[8].re);
    try std.testing.expectEqual(346, y5[8].im);

    const x6 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x6);
    const y6 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y6);

    @memcpy(x6.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y6.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .RowMajor, .Lower, n, k, alpha, A.ptr, k + 1, x6.ptr, -2, beta, y6.ptr, -2);

    try std.testing.expectEqual(16, y6[8].re);
    try std.testing.expectEqual(36, y6[8].im);
    try std.testing.expectEqual(34, y6[6].re);
    try std.testing.expectEqual(106, y6[6].im);
    try std.testing.expectEqual(52, y6[4].re);
    try std.testing.expectEqual(226, y6[4].im);
    try std.testing.expectEqual(70, y6[2].re);
    try std.testing.expectEqual(394, y6[2].im);
    try std.testing.expectEqual(-154, y6[0].re);
    try std.testing.expectEqual(346, y6[0].im);

    const x7 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x7);
    const y7 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y7);

    @memcpy(x7.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });
    @memcpy(y7.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .ColumnMajor, .Lower, n, k, alpha, A.ptr, k + 1, x7.ptr, 2, beta, y7.ptr, 2);

    try std.testing.expectEqual(11, y7[0].re);
    try std.testing.expectEqual(25, y7[0].im);
    try std.testing.expectEqual(29, y7[2].re);
    try std.testing.expectEqual(85, y7[2].im);
    try std.testing.expectEqual(47, y7[4].re);
    try std.testing.expectEqual(193, y7[4].im);
    try std.testing.expectEqual(65, y7[6].re);
    try std.testing.expectEqual(349, y7[6].im);
    try std.testing.expectEqual(-137, y7[8].re);
    try std.testing.expectEqual(313, y7[8].im);

    const x8 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x8);
    const y8 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y8);

    @memcpy(x8.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y8.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.hbmv(Complex(f64), .ColumnMajor, .Lower, n, k, alpha, A.ptr, k + 1, x8.ptr, -2, beta, y8.ptr, -2);

    try std.testing.expectEqual(11, y8[8].re);
    try std.testing.expectEqual(25, y8[8].im);
    try std.testing.expectEqual(29, y8[6].re);
    try std.testing.expectEqual(85, y8[6].im);
    try std.testing.expectEqual(47, y8[4].re);
    try std.testing.expectEqual(193, y8[4].im);
    try std.testing.expectEqual(65, y8[2].re);
    try std.testing.expectEqual(349, y8[2].im);
    try std.testing.expectEqual(-137, y8[0].re);
    try std.testing.expectEqual(313, y8[0].im);
}
