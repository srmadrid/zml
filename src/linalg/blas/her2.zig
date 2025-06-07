const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Scalar = types.Scalar;

pub inline fn her2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: Scalar(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (numericType) {
        .bool => @compileError("blas.her2 does not support bool."),
        .int, .float => @compileError("blas.her2 does not support int or float."),
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);
                    const t1 = T.init(alpha.re * y[@intCast(jy)].re + alpha.im * y[@intCast(jy)].im, -alpha.im * y[@intCast(jy)].re + alpha.re * y[@intCast(jy)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        A[@intCast(iaij)].re += t0.re * y[@intCast(iy)].re + t0.im * y[@intCast(iy)].im + t1.re * x[@intCast(ix)].re + t1.im * x[@intCast(ix)].im;
                        A[@intCast(iaij)].im += (t0.im * y[@intCast(iy)].re - t0.re * y[@intCast(iy)].im + t1.im * x[@intCast(ix)].re - t1.re * x[@intCast(ix)].im) * (-conj);

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    A[@intCast(iaij)].re += 2 * (t0.re * y[@intCast(jy)].re + t0.im * y[@intCast(jy)].im);
                    A[@intCast(iaij)].im = 0;

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
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);
                    const t1 = T.init(alpha.re * y[@intCast(jy)].re + alpha.im * y[@intCast(jy)].im, -alpha.im * y[@intCast(jy)].re + alpha.re * y[@intCast(jy)].im);
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)].re += 2 * (t0.re * y[@intCast(jy)].re + t0.im * y[@intCast(jy)].im);
                    A[@intCast(jaj)].im = 0;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        A[@intCast(iaij)].re += t0.re * y[@intCast(iy)].re + t0.im * y[@intCast(iy)].im + t1.re * x[@intCast(ix)].re + t1.im * x[@intCast(ix)].im;
                        A[@intCast(iaij)].im += (t0.im * y[@intCast(iy)].re - t0.re * y[@intCast(iy)].im + t1.im * x[@intCast(ix)].re - t1.re * x[@intCast(ix)].im) * (-conj);

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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.her2 only supports simple types."),
    }
}
