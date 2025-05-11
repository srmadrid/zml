const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Scalar = types.Scalar;

pub inline fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or k < 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: Scalar(T) = 1;
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
