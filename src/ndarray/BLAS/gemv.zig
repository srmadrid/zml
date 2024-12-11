const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Transpose = @import("../ndarray.zig").Transpose;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn gemv(comptime T: type, order: Order, trans: Transpose, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    var TRANS = trans;
    if (order == .RowMajor) {
        M = n;
        N = m;
        TRANS = if (trans == .NoTrans) .Trans else if (trans == .ConjNoTrans) .ConjTrans else if (trans == .Trans) .NoTrans else .ConjNoTrans;
    }

    if (lda < @max(1, M)) return;

    var LENX = N;
    var LENY = M;
    if (TRANS == .Trans or TRANS == .ConjTrans) {
        LENX = M;
        LENY = N;
    }

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.gemv does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
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

            if (TRANS == .NoTrans or TRANS == .ConjNoTrans) {
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
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    var t0: T = 0;

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
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

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            }
        },
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

            if (TRANS == .NoTrans) {
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
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                        y[@intCast(iy)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else if (TRANS == .ConjNoTrans) {
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
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.re * x[@intCast(jx)].im + alpha.im * x[@intCast(jx)].re);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        y[@intCast(iy)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                        y[@intCast(iy)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else if (TRANS == .Trans) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    var t0: T = T.init(0, 0);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
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

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    var t0: T = T.init(0, 0);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
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

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.gbmv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "gemv" {
    // TODO; Replace gbmv with gemv and make proper tests
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const kl = 1;
    const ku = 2;
    const alpha: f64 = 1;
    const beta: f64 = 1;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, m);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]f64{ 0, 0, 3, 7, 11, 0, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12, 0, 0 });
    @memcpy(x1.ptr, &[_]f64{ 1, 2, 3, 4, 5 });
    @memcpy(y1.ptr, &[_]f64{ 1, 2, 3, 4 });

    BLAS.gbmv(f64, .RowMajor, .NoTrans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x1.ptr, 1, beta, y1.ptr, 1);

    try std.testing.expectEqual(28, y1[0]);
    try std.testing.expectEqual(43, y1[1]);
    try std.testing.expectEqual(94, y1[2]);
    try std.testing.expectEqual(83, y1[3]);

    const x2 = try a.alloc(f64, n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * m);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]f64{ 5, 4, 3, 2, 1 });
    @memcpy(y2.ptr, &[_]f64{ 1, 0, 2, 0, 3, 0, 4, 0 });

    BLAS.gbmv(f64, .ColumnMajor, .NoTrans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x2.ptr, -1, beta, y2.ptr, 2);

    try std.testing.expectEqual(34, y2[0]);
    try std.testing.expectEqual(91, y2[2]);
    try std.testing.expectEqual(110, y2[4]);
    try std.testing.expectEqual(79, y2[6]);

    const x3 = try a.alloc(f64, m);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]f64{ 1, 2, 3, 4 });
    @memcpy(y3.ptr, &[_]f64{ 5, 0, 4, 0, 3, 0, 2, 0, 1, 0 });

    BLAS.gbmv(f64, .RowMajor, .Trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x3.ptr, 1, beta, y3.ptr, -2);

    try std.testing.expectEqual(23, y3[8]);
    try std.testing.expectEqual(35, y3[6]);
    try std.testing.expectEqual(92, y3[4]);
    try std.testing.expectEqual(71, y3[2]);
    try std.testing.expectEqual(20, y3[0]);

    const x4 = try a.alloc(f64, 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(f64, n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]f64{ 4, 0, 3, 0, 2, 0, 1, 0 });
    @memcpy(y4.ptr, &[_]f64{ 5, 4, 3, 2, 1 });

    BLAS.gbmv(f64, .ColumnMajor, .Trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x4.ptr, -2, beta, y4.ptr, -1);

    try std.testing.expectEqual(18, y4[4]);
    try std.testing.expectEqual(24, y4[3]);
    try std.testing.expectEqual(64, y4[2]);
    try std.testing.expectEqual(61, y4[1]);
    try std.testing.expectEqual(77, y4[0]);

    const alpha2 = Complex(f64).init(1, 1);
    const beta2 = Complex(f64).init(3, 3);

    const A2 = try a.alloc(Complex(f64), m * n);
    defer a.free(A2);
    const x5 = try a.alloc(Complex(f64), n);
    defer a.free(x5);
    const y5 = try a.alloc(Complex(f64), m);
    defer a.free(y5);

    @memcpy(A2.ptr, &[_]Complex(f64){
        Complex(f64).init(0, 0),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(7, 7),
        Complex(f64).init(11, 11),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(6, 6),
        Complex(f64).init(10, 10),
        Complex(f64).init(14, 14),
        Complex(f64).init(1, 1),
        Complex(f64).init(5, 5),
        Complex(f64).init(9, 9),
        Complex(f64).init(13, 13),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(8, 8),
        Complex(f64).init(12, 12),
        Complex(f64).init(0, 0),
        Complex(f64).init(0, 0),
    });
    @memcpy(x5.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });
    @memcpy(y5.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
    });

    BLAS.gbmv(Complex(f64), .RowMajor, .NoTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x5.ptr, 1, beta2, y5.ptr, 1);

    try std.testing.expectEqual(-111, y5[0].re);
    try std.testing.expectEqual(97, y5[0].im);
    try std.testing.expectEqual(-164, y5[1].re);
    try std.testing.expectEqual(144, y5[1].im);
    try std.testing.expectEqual(-367, y5[2].re);
    try std.testing.expectEqual(313, y5[2].im);
    try std.testing.expectEqual(-316, y5[3].re);
    try std.testing.expectEqual(290, y5[3].im);

    const x6 = try a.alloc(Complex(f64), n);
    defer a.free(x6);
    const y6 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y6);

    @memcpy(x6.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });
    @memcpy(y6.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
    });

    BLAS.gbmv(Complex(f64), .ColumnMajor, .NoTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x6.ptr, -1, beta2, y6.ptr, 2);

    try std.testing.expectEqual(-135, y6[0].re);
    try std.testing.expectEqual(115, y6[0].im);
    try std.testing.expectEqual(-356, y6[2].re);
    try std.testing.expectEqual(310, y6[2].im);
    try std.testing.expectEqual(-431, y6[4].re);
    try std.testing.expectEqual(381, y6[4].im);
    try std.testing.expectEqual(-300, y6[6].re);
    try std.testing.expectEqual(284, y6[6].im);

    const x7 = try a.alloc(Complex(f64), m);
    defer a.free(x7);
    const y7 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y7);

    @memcpy(x7.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
    });
    @memcpy(y7.ptr, &[_]Complex(f64){
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

    BLAS.gbmv(Complex(f64), .RowMajor, .Trans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x7.ptr, 1, beta2, y7.ptr, -2);

    try std.testing.expectEqual(-91, y7[8].re);
    try std.testing.expectEqual(75, y7[8].im);
    try std.testing.expectEqual(-132, y7[6].re);
    try std.testing.expectEqual(124, y7[6].im);
    try std.testing.expectEqual(-359, y7[4].re);
    try std.testing.expectEqual(301, y7[4].im);
    try std.testing.expectEqual(-268, y7[2].re);
    try std.testing.expectEqual(246, y7[2].im);
    try std.testing.expectEqual(-63, y7[0].re);
    try std.testing.expectEqual(59, y7[0].im);

    const x8 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x8);
    const y8 = try a.alloc(Complex(f64), n);
    defer a.free(y8);

    @memcpy(x8.ptr, &[_]Complex(f64){
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
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.gbmv(Complex(f64), .ColumnMajor, .Trans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x8.ptr, -2, beta2, y8.ptr, -1);

    try std.testing.expectEqual(-71, y8[4].re);
    try std.testing.expectEqual(57, y8[4].im);
    try std.testing.expectEqual(-88, y8[3].re);
    try std.testing.expectEqual(90, y8[3].im);
    try std.testing.expectEqual(-247, y8[2].re);
    try std.testing.expectEqual(193, y8[2].im);
    try std.testing.expectEqual(-228, y8[1].re);
    try std.testing.expectEqual(202, y8[1].im);
    try std.testing.expectEqual(-291, y8[0].re);
    try std.testing.expectEqual(257, y8[0].im);

    const x9 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x9);
    const y9 = try a.alloc(Complex(f64), m);
    defer a.free(y9);

    @memcpy(x9.ptr, &[_]Complex(f64){
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
    @memcpy(y9.ptr, &[_]Complex(f64){
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.gbmv(Complex(f64), .RowMajor, .ConjNoTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x9.ptr, -2, beta2, y9.ptr, -1);

    try std.testing.expectEqual(85, y9[3].re);
    try std.testing.expectEqual(117, y9[3].im);
    try std.testing.expectEqual(126, y9[2].re);
    try std.testing.expectEqual(182, y9[2].im);
    try std.testing.expectEqual(301, y9[1].re);
    try std.testing.expectEqual(373, y9[1].im);
    try std.testing.expectEqual(272, y9[0].re);
    try std.testing.expectEqual(334, y9[0].im);

    const x10 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x10);
    const y10 = try a.alloc(Complex(f64), m);
    defer a.free(y10);

    @memcpy(x10.ptr, &[_]Complex(f64){
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
    @memcpy(y10.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
    });

    BLAS.gbmv(Complex(f64), .ColumnMajor, .ConjNoTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x10.ptr, 2, beta2, y10.ptr, 1);

    try std.testing.expectEqual(103, y10[0].re);
    try std.testing.expectEqual(141, y10[0].im);
    try std.testing.expectEqual(292, y10[1].re);
    try std.testing.expectEqual(374, y10[1].im);
    try std.testing.expectEqual(369, y10[2].re);
    try std.testing.expectEqual(437, y10[2].im);
    try std.testing.expectEqual(266, y10[3].re);
    try std.testing.expectEqual(318, y10[3].im);

    const x11 = try a.alloc(Complex(f64), m);
    defer a.free(x11);
    const y11 = try a.alloc(Complex(f64), n);
    defer a.free(y11);

    @memcpy(x11.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });
    @memcpy(y11.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.gbmv(Complex(f64), .RowMajor, .ConjTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x11.ptr, -1, beta2, y11.ptr, -1);

    try std.testing.expectEqual(63, y11[4].re);
    try std.testing.expectEqual(97, y11[4].im);
    try std.testing.expectEqual(106, y11[3].re);
    try std.testing.expectEqual(150, y11[3].im);
    try std.testing.expectEqual(289, y11[2].re);
    try std.testing.expectEqual(365, y11[2].im);
    try std.testing.expectEqual(228, y11[1].re);
    try std.testing.expectEqual(286, y11[1].im);
    try std.testing.expectEqual(47, y11[0].re);
    try std.testing.expectEqual(69, y11[0].im);

    const x12 = try a.alloc(Complex(f64), m);
    defer a.free(x12);
    const y12 = try a.alloc(Complex(f64), n);
    defer a.free(y12);

    @memcpy(x12.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
    });
    @memcpy(y12.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.gbmv(Complex(f64), .ColumnMajor, .ConjTrans, m, n, kl, ku, alpha2, A2.ptr, kl + ku + 1, x12.ptr, 1, beta2, y12.ptr, 1);

    try std.testing.expectEqual(45, y12[0].re);
    try std.testing.expectEqual(77, y12[0].im);
    try std.testing.expectEqual(72, y12[1].re);
    try std.testing.expectEqual(106, y12[1].im);
    try std.testing.expectEqual(181, y12[2].re);
    try std.testing.expectEqual(253, y12[2].im);
    try std.testing.expectEqual(184, y12[3].re);
    try std.testing.expectEqual(246, y12[3].im);
    try std.testing.expectEqual(245, y12[4].re);
    try std.testing.expectEqual(297, y12[4].im);
}
