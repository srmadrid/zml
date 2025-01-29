const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or k < 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: scalar(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < k + 1) return;

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.hbmv does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("BLAS.hbmv does not support int or float."),
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
                } else if (beta.re != 1 and beta.im != 0) {
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
                } else if (beta.re != 1 and beta.im != 0) {
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
                } else if (beta.re != 1 and beta.im != 0) {
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
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.hbmv only supports simple types."),
        .Unsupported => unreachable,
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
    const x1 = try a.alloc(Complex(f64), n);
    defer a.free(x1);
    const y1 = try a.alloc(Complex(f64), n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]Complex(f64){
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 3),
        Complex(f64).init(4, -1),
        Complex(f64).init(5, 2),
        Complex(f64).init(6, -3),
        Complex(f64).init(1, 0),
        Complex(f64).init(7, 0),
        Complex(f64).init(11, 0),
        Complex(f64).init(14, 0),
        Complex(f64).init(16, 0),
        Complex(f64).init(2, -3),
        Complex(f64).init(4, 1),
        Complex(f64).init(5, -2),
        Complex(f64).init(6, 3),
        Complex(f64).init(0, 0),
    });
    @memcpy(x1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });
    @memcpy(y1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.hbmv(Complex(f64), .RowMajor, .Upper, n, k, alpha, A.ptr, k + 1, x1.ptr, 1, beta, y1.ptr, 1);

    try std.testing.expectEqual(-26, y1[0].re);
    try std.testing.expectEqual(20, y1[0].im);
    try std.testing.expectEqual(-24, y1[1].re);
    try std.testing.expectEqual(108, y1[1].im);
    try std.testing.expectEqual(-1, y1[2].re);
    try std.testing.expectEqual(127, y1[2].im);
    try std.testing.expectEqual(-19, y1[3].re);
    try std.testing.expectEqual(343, y1[3].im);
    try std.testing.expectEqual(-28, y1[4].re);
    try std.testing.expectEqual(440, y1[4].im);

    const x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);
    const y2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });
    @memcpy(y2.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });

    BLAS.hbmv(Complex(f64), .ColumnMajor, .Upper, n, k, alpha, A.ptr, k + 1, x2.ptr, -1, beta, y2.ptr, 2);

    try std.testing.expectEqual(-2, y2[0].re);
    try std.testing.expectEqual(44, y2[0].im);
    try std.testing.expectEqual(15, y2[2].re);
    try std.testing.expectEqual(133, y2[2].im);
    try std.testing.expectEqual(-38, y2[4].re);
    try std.testing.expectEqual(164, y2[4].im);
    try std.testing.expectEqual(-32, y2[6].re);
    try std.testing.expectEqual(526, y2[6].im);
    try std.testing.expectEqual(-33, y2[8].re);
    try std.testing.expectEqual(523, y2[8].im);

    const x3 = try a.alloc(Complex(f64), n);
    defer a.free(x3);
    const y3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });
    @memcpy(y3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });

    BLAS.hbmv(Complex(f64), .RowMajor, .Lower, n, k, alpha, A.ptr, k + 1, x3.ptr, 1, beta, y3.ptr, -2);

    try std.testing.expectEqual(-16, y3[8].re);
    try std.testing.expectEqual(42, y3[8].im);
    try std.testing.expectEqual(-45, y3[6].re);
    try std.testing.expectEqual(129, y3[6].im);
    try std.testing.expectEqual(4, y3[4].re);
    try std.testing.expectEqual(170, y3[4].im);
    try std.testing.expectEqual(-32, y3[2].re);
    try std.testing.expectEqual(526, y3[2].im);
    try std.testing.expectEqual(-33, y3[0].re);
    try std.testing.expectEqual(523, y3[0].im);

    const x4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(Complex(f64), n);
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
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.hbmv(Complex(f64), .ColumnMajor, .Lower, n, k, alpha, A.ptr, k + 1, x4.ptr, -2, beta, y4.ptr, -1);

    try std.testing.expectEqual(16, y4[4].re);
    try std.testing.expectEqual(26, y4[4].im);
    try std.testing.expectEqual(2, y4[3].re);
    try std.testing.expectEqual(106, y4[3].im);
    try std.testing.expectEqual(-29, y4[2].re);
    try std.testing.expectEqual(123, y4[2].im);
    try std.testing.expectEqual(-19, y4[1].re);
    try std.testing.expectEqual(343, y4[1].im);
    try std.testing.expectEqual(-28, y4[0].re);
    try std.testing.expectEqual(440, y4[0].im);
}
