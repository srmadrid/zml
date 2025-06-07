const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

pub inline fn symv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (numericType) {
        .bool => @compileError("blas.symv does not support bool."),
        .int, .float => {
            if (alpha == 0 and beta == 1) return;

            if (alpha == 0) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                return;
            }

            if (UPLO == .Upper) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    var t1: T = 0;

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        t1 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)] += t0 * A[@intCast(iaij)] + alpha * t1;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            } else {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    var t1: T = 0;
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        t1 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)] += t0 * A[@intCast(jaj)] + alpha * t1;

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .cfloat => @compileError("blas.symv does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.symv only supports simple types."),
    }
}
