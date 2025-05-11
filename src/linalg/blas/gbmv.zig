const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Transpose = blas.Transpose;

pub inline fn gbmv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (m <= 0 or n <= 0 or kl < 0 or ku < 0 or lda < kl + ku + 1) return;

    var M = m;
    var N = n;
    var KL = kl;
    var KU = ku;
    var TRANSA = transA;
    if (order == .RowMajor) {
        M = n;
        N = m;
        KL = ku;
        KU = kl;
        TRANSA = if (transA == .NoTrans) .Trans else if (transA == .ConjNoTrans) .ConjTrans else if (transA == .Trans) .NoTrans else .ConjNoTrans;
    }

    var LENX = N;
    var LENY = M;
    if (TRANSA == .Trans or TRANSA == .ConjTrans) {
        LENX = M;
        LENY = N;
    }

    switch (numericType) {
        .bool => @compileError("blas.gbmv does not support bool."),
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

            if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
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
                var ky: isize = 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy + ky else ky;
                    while (i <= I1) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    if (j >= KU) ky += incy;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                var kx: isize = 0;
                while (j < N) {
                    var t0: T = 0;
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                    while (i <= I1) {
                        t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    if (beta == 0) {
                        y[@intCast(jy)] = 0;
                    } else if (beta != 1) {
                        y[@intCast(jy)] *= beta;
                    }

                    y[@intCast(jy)] += alpha * t0;

                    if (j >= KU) kx += incx;

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            }
        },
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

            if (TRANSA == .NoTrans) {
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
                var ky: isize = 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy + ky else ky;
                    while (i <= I1) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                        y[@intCast(iy)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    if (j >= KU) ky += incy;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else if (TRANSA == .ConjNoTrans) {
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
                var ky: isize = 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy + ky else ky;
                    while (i <= I1) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                        y[@intCast(iy)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    if (j >= KU) ky += incy;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else if (TRANSA == .Trans) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                var kx: isize = 0;
                while (j < N) {
                    var t0: T = T.init(0, 0);
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                    while (i <= I1) {
                        t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                        t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    if (beta.re == 0 and beta.im == 0) {
                        y[@intCast(jy)].re = 0;
                        y[@intCast(jy)].im = 0;
                    } else if (beta.re != 1 or beta.im != 0) {
                        const tmp = y[@intCast(jy)].re * beta.re - y[@intCast(jy)].im * beta.im;
                        y[@intCast(jy)].im = y[@intCast(jy)].re * beta.im + y[@intCast(jy)].im * beta.re;
                        y[@intCast(jy)].re = tmp;
                    }

                    y[@intCast(jy)].re += alpha.re * t0.re - alpha.im * t0.im;
                    y[@intCast(jy)].im += alpha.re * t0.im + alpha.im * t0.re;

                    if (j >= KU) kx += incx;

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                var kx: isize = 0;
                while (j < N) {
                    var t0: T = T.init(0, 0);
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                    while (i <= I1) {
                        t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                        t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    if (beta.re == 0 and beta.im == 0) {
                        y[@intCast(jy)].re = 0;
                        y[@intCast(jy)].im = 0;
                    } else if (beta.re != 1 or beta.im != 0) {
                        const tmp = y[@intCast(jy)].re * beta.re - y[@intCast(jy)].im * beta.im;
                        y[@intCast(jy)].im = y[@intCast(jy)].re * beta.im + y[@intCast(jy)].im * beta.re;
                        y[@intCast(jy)].re = tmp;
                    }

                    y[@intCast(jy)].re += alpha.re * t0.re - alpha.im * t0.im;
                    y[@intCast(jy)].im += alpha.re * t0.im + alpha.im * t0.re;

                    if (j >= KU) kx += incx;

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.gbmv only supports simple types."),
        .unsupported => unreachable,
    }
}
