const std = @import("std");
const core = @import("../../core/core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const scalar = core.supported.scalar;

pub inline fn hemv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: scalar(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("blas.hemv does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("blas.hemv does not support int or float."),
        .Complex => {
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
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    var t1 = T.init(0, 0);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im * conj;
                        y[@intCast(iy)].im += t0.im * A[@intCast(iaij)].re + t0.re * A[@intCast(iaij)].im * conj;

                        t1.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * conj * x[@intCast(ix)].im;
                        t1.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * conj * x[@intCast(ix)].re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)].re += t0.re * A[@intCast(iaij)].re + alpha.re * t1.re - alpha.im * t1.im;
                    y[@intCast(jy)].im += t0.im * A[@intCast(iaij)].re + alpha.re * t1.im + alpha.im * t1.re;

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

                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    var t1 = T.init(0, 0);
                    const I0: isize = j + 1;
                    const I1: isize = N;

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
                    jaj += ldap1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.hemv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "hemv" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const alpha = Complex(f64).init(1, 1);
    const beta = Complex(f64).init(3, 3);

    const A = try a.alloc(Complex(f64), n * n);
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
        Complex(f64).init(16, 16),
        Complex(f64).init(17, 17),
        Complex(f64).init(18, 18),
        Complex(f64).init(19, 19),
        Complex(f64).init(20, 20),
        Complex(f64).init(21, 21),
        Complex(f64).init(22, 22),
        Complex(f64).init(23, 23),
        Complex(f64).init(24, 24),
        Complex(f64).init(25, 25),
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

    blas.hemv(Complex(f64), .RowMajor, .Upper, n, alpha, A.ptr, n, x1.ptr, 2, beta, y1.ptr, 2);

    try std.testing.expectEqual(-217, y1[0].re);
    try std.testing.expectEqual(197, y1[0].im);
    try std.testing.expectEqual(-443, y1[2].re);
    try std.testing.expectEqual(455, y1[2].im);
    try std.testing.expectEqual(-483, y1[4].re);
    try std.testing.expectEqual(703, y1[4].im);
    try std.testing.expectEqual(-217, y1[6].re);
    try std.testing.expectEqual(925, y1[6].im);
    try std.testing.expectEqual(475, y1[8].re);
    try std.testing.expectEqual(1105, y1[8].im);

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

    blas.hemv(Complex(f64), .RowMajor, .Upper, n, alpha, A.ptr, n, x2.ptr, -2, beta, y2.ptr, -2);

    try std.testing.expectEqual(-217, y2[8].re);
    try std.testing.expectEqual(197, y2[8].im);
    try std.testing.expectEqual(-443, y2[6].re);
    try std.testing.expectEqual(455, y2[6].im);
    try std.testing.expectEqual(-483, y2[4].re);
    try std.testing.expectEqual(703, y2[4].im);
    try std.testing.expectEqual(-217, y2[2].re);
    try std.testing.expectEqual(925, y2[2].im);
    try std.testing.expectEqual(475, y2[0].re);
    try std.testing.expectEqual(1105, y2[0].im);

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

    blas.hemv(Complex(f64), .ColumnMajor, .Upper, n, alpha, A.ptr, n, x3.ptr, 2, beta, y3.ptr, 2);

    try std.testing.expectEqual(-857, y3[0].re);
    try std.testing.expectEqual(757, y3[0].im);
    try std.testing.expectEqual(-851, y3[2].re);
    try std.testing.expectEqual(839, y3[2].im);
    try std.testing.expectEqual(-667, y3[4].re);
    try std.testing.expectEqual(967, y3[4].im);
    try std.testing.expectEqual(-185, y3[6].re);
    try std.testing.expectEqual(1157, y3[6].im);
    try std.testing.expectEqual(715, y3[8].re);
    try std.testing.expectEqual(1425, y3[8].im);

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

    blas.hemv(Complex(f64), .ColumnMajor, .Upper, n, alpha, A.ptr, n, x4.ptr, -2, beta, y4.ptr, -2);

    try std.testing.expectEqual(-857, y4[8].re);
    try std.testing.expectEqual(757, y4[8].im);
    try std.testing.expectEqual(-851, y4[6].re);
    try std.testing.expectEqual(839, y4[6].im);
    try std.testing.expectEqual(-667, y4[4].re);
    try std.testing.expectEqual(967, y4[4].im);
    try std.testing.expectEqual(-185, y4[2].re);
    try std.testing.expectEqual(1157, y4[2].im);
    try std.testing.expectEqual(715, y4[0].re);
    try std.testing.expectEqual(1425, y4[0].im);

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

    blas.hemv(Complex(f64), .RowMajor, .Lower, n, alpha, A.ptr, n, x5.ptr, 2, beta, y5.ptr, 2);

    try std.testing.expectEqual(747, y5[0].re);
    try std.testing.expectEqual(865, y5[0].im);
    try std.testing.expectEqual(723, y5[2].re);
    try std.testing.expectEqual(929, y5[2].im);
    try std.testing.expectEqual(513, y5[4].re);
    try std.testing.expectEqual(1003, y5[4].im);
    try std.testing.expectEqual(-3, y5[6].re);
    try std.testing.expectEqual(1103, y5[6].im);
    try std.testing.expectEqual(-945, y5[8].re);
    try std.testing.expectEqual(1245, y5[8].im);

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

    blas.hemv(Complex(f64), .RowMajor, .Lower, n, alpha, A.ptr, n, x6.ptr, -2, beta, y6.ptr, -2);

    try std.testing.expectEqual(747, y6[8].re);
    try std.testing.expectEqual(865, y6[8].im);
    try std.testing.expectEqual(723, y6[6].re);
    try std.testing.expectEqual(929, y6[6].im);
    try std.testing.expectEqual(513, y6[4].re);
    try std.testing.expectEqual(1003, y6[4].im);
    try std.testing.expectEqual(-3, y6[2].re);
    try std.testing.expectEqual(1103, y6[2].im);
    try std.testing.expectEqual(-945, y6[0].re);
    try std.testing.expectEqual(1245, y6[0].im);

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

    blas.hemv(Complex(f64), .ColumnMajor, .Lower, n, alpha, A.ptr, n, x7.ptr, 2, beta, y7.ptr, 2);

    try std.testing.expectEqual(187, y7[0].re);
    try std.testing.expectEqual(225, y7[0].im);
    try std.testing.expectEqual(371, y7[2].re);
    try std.testing.expectEqual(505, y7[2].im);
    try std.testing.expectEqual(377, y7[4].re);
    try std.testing.expectEqual(739, y7[4].im);
    try std.testing.expectEqual(85, y7[6].re);
    try std.testing.expectEqual(911, y7[6].im);
    try std.testing.expectEqual(-625, y7[8].re);
    try std.testing.expectEqual(1005, y7[8].im);

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

    blas.hemv(Complex(f64), .ColumnMajor, .Lower, n, alpha, A.ptr, n, x8.ptr, -2, beta, y8.ptr, -2);

    try std.testing.expectEqual(187, y8[8].re);
    try std.testing.expectEqual(225, y8[8].im);
    try std.testing.expectEqual(371, y8[6].re);
    try std.testing.expectEqual(505, y8[6].im);
    try std.testing.expectEqual(377, y8[4].re);
    try std.testing.expectEqual(739, y8[4].im);
    try std.testing.expectEqual(85, y8[2].re);
    try std.testing.expectEqual(911, y8[2].im);
    try std.testing.expectEqual(-625, y8[0].re);
    try std.testing.expectEqual(1005, y8[0].im);
}
