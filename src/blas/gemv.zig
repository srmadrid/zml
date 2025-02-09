const std = @import("std");
const core = @import("../core/core.zig");
const Order = @import("../ndarray/ndarray.zig").Order;
const Transpose = @import("../ndarray/ndarray.zig").Transpose;
const blas = @import("blas.zig");

const scalar = core.supported.scalar;

pub inline fn gemv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    var TRANSA = transA;
    if (order == .RowMajor) {
        M = n;
        N = m;
        TRANSA = if (transA == .NoTrans) .Trans else if (transA == .ConjNoTrans) .ConjTrans else if (transA == .Trans) .NoTrans else .ConjNoTrans;
    }

    if (lda < @max(1, M)) return;

    var LENX = N;
    var LENY = M;
    if (TRANSA == .Trans or TRANSA == .ConjTrans) {
        LENX = M;
        LENY = N;
    }

    switch (supported) {
        .BuiltinBool => @compileError("blas.gemv does not support bool."),
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

            if (TRANSA == .NoTrans) {
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
            } else if (TRANSA == .ConjNoTrans) {
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
            } else if (TRANSA == .Trans) {
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
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.gemv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "gemv" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const alpha: f64 = 2;
    const beta: f64 = 3;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, 2 * m);
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
    });

    blas.gemv(f64, .RowMajor, .NoTrans, m, n, alpha, A.ptr, n, x1.ptr, 2, beta, y1.ptr, 2);

    try std.testing.expectEqual(113, y1[0]);
    try std.testing.expectEqual(266, y1[2]);
    try std.testing.expectEqual(419, y1[4]);
    try std.testing.expectEqual(572, y1[6]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * m);
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
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.gemv(f64, .RowMajor, .NoTrans, m, n, alpha, A.ptr, n, x2.ptr, -2, beta, y2.ptr, -2);

    try std.testing.expectEqual(113, y2[6]);
    try std.testing.expectEqual(266, y2[4]);
    try std.testing.expectEqual(419, y2[2]);
    try std.testing.expectEqual(572, y2[0]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * m);
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
    });

    blas.gemv(f64, .ColumnMajor, .NoTrans, m, n, alpha, A.ptr, m, x3.ptr, 2, beta, y3.ptr, 2);

    try std.testing.expectEqual(353, y3[0]);
    try std.testing.expectEqual(386, y3[2]);
    try std.testing.expectEqual(419, y3[4]);
    try std.testing.expectEqual(452, y3[6]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * m);
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
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    blas.gemv(f64, .ColumnMajor, .NoTrans, m, n, alpha, A.ptr, m, x4.ptr, -2, beta, y4.ptr, -2);

    try std.testing.expectEqual(353, y4[6]);
    try std.testing.expectEqual(386, y4[4]);
    try std.testing.expectEqual(419, y4[2]);
    try std.testing.expectEqual(452, y4[0]);

    const x5 = try a.alloc(f64, 2 * m);
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

    blas.gemv(f64, .RowMajor, .Trans, m, n, alpha, A.ptr, n, x5.ptr, 2, beta, y5.ptr, 2);

    try std.testing.expectEqual(223, y5[0]);
    try std.testing.expectEqual(246, y5[2]);
    try std.testing.expectEqual(269, y5[4]);
    try std.testing.expectEqual(292, y5[6]);
    try std.testing.expectEqual(315, y5[8]);

    const x6 = try a.alloc(f64, 2 * m);
    defer a.free(x6);
    const y6 = try a.alloc(f64, 2 * n);
    defer a.free(y6);

    @memcpy(x6.ptr, &[_]f64{
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

    blas.gemv(f64, .RowMajor, .Trans, m, n, alpha, A.ptr, n, x6.ptr, -2, beta, y6.ptr, -2);

    try std.testing.expectEqual(223, y6[8]);
    try std.testing.expectEqual(246, y6[6]);
    try std.testing.expectEqual(269, y6[4]);
    try std.testing.expectEqual(292, y6[2]);
    try std.testing.expectEqual(315, y6[0]);

    const x7 = try a.alloc(f64, 2 * m);
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

    blas.gemv(f64, .ColumnMajor, .Trans, m, n, alpha, A.ptr, m, x7.ptr, 2, beta, y7.ptr, 2);

    try std.testing.expectEqual(63, y7[0]);
    try std.testing.expectEqual(146, y7[2]);
    try std.testing.expectEqual(229, y7[4]);
    try std.testing.expectEqual(312, y7[6]);
    try std.testing.expectEqual(395, y7[8]);

    const x8 = try a.alloc(f64, 2 * m);
    defer a.free(x8);
    const y8 = try a.alloc(f64, 2 * n);
    defer a.free(y8);

    @memcpy(x8.ptr, &[_]f64{
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

    blas.gemv(f64, .ColumnMajor, .Trans, m, n, alpha, A.ptr, m, x8.ptr, -2, beta, y8.ptr, -2);

    try std.testing.expectEqual(63, y8[8]);
    try std.testing.expectEqual(146, y8[6]);
    try std.testing.expectEqual(229, y8[4]);
    try std.testing.expectEqual(312, y8[2]);
    try std.testing.expectEqual(395, y8[0]);

    const gamma = Complex(f64).init(2, 2);
    const delta = Complex(f64).init(3, 3);

    const B = try a.alloc(Complex(f64), m * n);
    defer a.free(B);
    const x9 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x9);
    const y9 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y9);

    @memcpy(B.ptr, &[_]Complex(f64){
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
    });
    @memcpy(x9.ptr, &[_]Complex(f64){
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
    @memcpy(y9.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .RowMajor, .NoTrans, m, n, gamma, B.ptr, n, x9.ptr, 2, delta, y9.ptr, 2);

    try std.testing.expectEqual(-440, y9[0].re);
    try std.testing.expectEqual(386, y9[0].im);
    try std.testing.expectEqual(-1040, y9[2].re);
    try std.testing.expectEqual(892, y9[2].im);
    try std.testing.expectEqual(-1640, y9[4].re);
    try std.testing.expectEqual(1398, y9[4].im);
    try std.testing.expectEqual(-2240, y9[6].re);
    try std.testing.expectEqual(1904, y9[6].im);

    const x10 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x10);
    const y10 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y10);

    @memcpy(x10.ptr, &[_]Complex(f64){
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
    @memcpy(y10.ptr, &[_]Complex(f64){
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .RowMajor, .NoTrans, m, n, gamma, B.ptr, n, x10.ptr, -2, delta, y10.ptr, -2);

    try std.testing.expectEqual(-440, y10[6].re);
    try std.testing.expectEqual(386, y10[6].im);
    try std.testing.expectEqual(-1040, y10[4].re);
    try std.testing.expectEqual(892, y10[4].im);
    try std.testing.expectEqual(-1640, y10[2].re);
    try std.testing.expectEqual(1398, y10[2].im);
    try std.testing.expectEqual(-2240, y10[0].re);
    try std.testing.expectEqual(1904, y10[0].im);

    const x11 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x11);
    const y11 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y11);

    @memcpy(x11.ptr, &[_]Complex(f64){
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
    @memcpy(y11.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .ColumnMajor, .NoTrans, m, n, gamma, B.ptr, m, x11.ptr, 2, delta, y11.ptr, 2);

    try std.testing.expectEqual(-1400, y11[0].re);
    try std.testing.expectEqual(1226, y11[0].im);
    try std.testing.expectEqual(-1520, y11[2].re);
    try std.testing.expectEqual(1332, y11[2].im);
    try std.testing.expectEqual(-1640, y11[4].re);
    try std.testing.expectEqual(1438, y11[4].im);
    try std.testing.expectEqual(-1760, y11[6].re);
    try std.testing.expectEqual(1544, y11[6].im);

    const x12 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x12);
    const y12 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y12);

    @memcpy(x12.ptr, &[_]Complex(f64){
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
    @memcpy(y12.ptr, &[_]Complex(f64){
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .ColumnMajor, .NoTrans, m, n, gamma, B.ptr, m, x12.ptr, -2, delta, y12.ptr, -2);

    try std.testing.expectEqual(-1400, y12[6].re);
    try std.testing.expectEqual(1226, y12[6].im);
    try std.testing.expectEqual(-1520, y12[4].re);
    try std.testing.expectEqual(1332, y12[4].im);
    try std.testing.expectEqual(-1640, y12[2].re);
    try std.testing.expectEqual(1438, y12[2].im);
    try std.testing.expectEqual(-1760, y12[0].re);
    try std.testing.expectEqual(1544, y12[0].im);

    const x13 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x13);
    const y13 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y13);

    @memcpy(x13.ptr, &[_]Complex(f64){
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
    @memcpy(y13.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .RowMajor, .ConjNoTrans, m, n, gamma, B.ptr, n, x13.ptr, 2, delta, y13.ptr, 2);

    try std.testing.expectEqual(380, y13[0].re);
    try std.testing.expectEqual(446, y13[0].im);
    try std.testing.expectEqual(880, y13[2].re);
    try std.testing.expectEqual(1052, y13[2].im);
    try std.testing.expectEqual(1380, y13[4].re);
    try std.testing.expectEqual(1658, y13[4].im);
    try std.testing.expectEqual(1880, y13[6].re);
    try std.testing.expectEqual(2264, y13[6].im);

    const x14 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x14);
    const y14 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y14);

    @memcpy(x14.ptr, &[_]Complex(f64){
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
    @memcpy(y14.ptr, &[_]Complex(f64){
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .RowMajor, .ConjNoTrans, m, n, gamma, B.ptr, n, x14.ptr, -2, delta, y14.ptr, -2);

    try std.testing.expectEqual(380, y14[6].re);
    try std.testing.expectEqual(446, y14[6].im);
    try std.testing.expectEqual(880, y14[4].re);
    try std.testing.expectEqual(1052, y14[4].im);
    try std.testing.expectEqual(1380, y14[2].re);
    try std.testing.expectEqual(1658, y14[2].im);
    try std.testing.expectEqual(1880, y14[0].re);
    try std.testing.expectEqual(2264, y14[0].im);

    const x15 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x15);
    const y15 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y15);

    @memcpy(x15.ptr, &[_]Complex(f64){
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
    @memcpy(y15.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .ColumnMajor, .ConjNoTrans, m, n, gamma, B.ptr, m, x15.ptr, 2, delta, y15.ptr, 2);

    try std.testing.expectEqual(1220, y15[0].re);
    try std.testing.expectEqual(1406, y15[0].im);
    try std.testing.expectEqual(1320, y15[2].re);
    try std.testing.expectEqual(1532, y15[2].im);
    try std.testing.expectEqual(1420, y15[4].re);
    try std.testing.expectEqual(1658, y15[4].im);
    try std.testing.expectEqual(1520, y15[6].re);
    try std.testing.expectEqual(1784, y15[6].im);

    const x16 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x16);
    const y16 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(y16);

    @memcpy(x16.ptr, &[_]Complex(f64){
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
    @memcpy(y16.ptr, &[_]Complex(f64){
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.gemv(Complex(f64), .ColumnMajor, .ConjNoTrans, m, n, gamma, B.ptr, m, x16.ptr, -2, delta, y16.ptr, -2);

    try std.testing.expectEqual(1220, y16[6].re);
    try std.testing.expectEqual(1406, y16[6].im);
    try std.testing.expectEqual(1320, y16[4].re);
    try std.testing.expectEqual(1532, y16[4].im);
    try std.testing.expectEqual(1420, y16[2].re);
    try std.testing.expectEqual(1658, y16[2].im);
    try std.testing.expectEqual(1520, y16[0].re);
    try std.testing.expectEqual(1784, y16[0].im);

    const x17 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x17);
    const y17 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y17);

    @memcpy(x17.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y17.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .RowMajor, .Trans, m, n, gamma, B.ptr, n, x17.ptr, 2, delta, y17.ptr, 2);

    try std.testing.expectEqual(-880, y17[0].re);
    try std.testing.expectEqual(750, y17[0].im);
    try std.testing.expectEqual(-960, y17[2].re);
    try std.testing.expectEqual(820, y17[2].im);
    try std.testing.expectEqual(-1040, y17[4].re);
    try std.testing.expectEqual(890, y17[4].im);
    try std.testing.expectEqual(-1120, y17[6].re);
    try std.testing.expectEqual(960, y17[6].im);
    try std.testing.expectEqual(-1200, y17[8].re);
    try std.testing.expectEqual(1030, y17[8].im);

    const x18 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x18);
    const y18 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y18);

    @memcpy(x18.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y18.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .RowMajor, .Trans, m, n, gamma, B.ptr, n, x18.ptr, -2, delta, y18.ptr, -2);

    try std.testing.expectEqual(-880, y18[8].re);
    try std.testing.expectEqual(750, y18[8].im);
    try std.testing.expectEqual(-960, y18[6].re);
    try std.testing.expectEqual(820, y18[6].im);
    try std.testing.expectEqual(-1040, y18[4].re);
    try std.testing.expectEqual(890, y18[4].im);
    try std.testing.expectEqual(-1120, y18[2].re);
    try std.testing.expectEqual(960, y18[2].im);
    try std.testing.expectEqual(-1200, y18[0].re);
    try std.testing.expectEqual(1030, y18[0].im);

    const x19 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x19);
    const y19 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y19);

    @memcpy(x19.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y19.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .ColumnMajor, .Trans, m, n, gamma, B.ptr, m, x19.ptr, 2, delta, y19.ptr, 2);

    try std.testing.expectEqual(-240, y19[0].re);
    try std.testing.expectEqual(206, y19[0].im);
    try std.testing.expectEqual(-560, y19[2].re);
    try std.testing.expectEqual(468, y19[2].im);
    try std.testing.expectEqual(-880, y19[4].re);
    try std.testing.expectEqual(730, y19[4].im);
    try std.testing.expectEqual(-1200, y19[6].re);
    try std.testing.expectEqual(992, y19[6].im);
    try std.testing.expectEqual(-1520, y19[8].re);
    try std.testing.expectEqual(1254, y19[8].im);

    const x20 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x20);
    const y20 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y20);

    @memcpy(x20.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y20.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .ColumnMajor, .Trans, m, n, gamma, B.ptr, m, x20.ptr, -2, delta, y20.ptr, -2);

    try std.testing.expectEqual(-240, y20[8].re);
    try std.testing.expectEqual(206, y20[8].im);
    try std.testing.expectEqual(-560, y20[6].re);
    try std.testing.expectEqual(468, y20[6].im);
    try std.testing.expectEqual(-880, y20[4].re);
    try std.testing.expectEqual(730, y20[4].im);
    try std.testing.expectEqual(-1200, y20[2].re);
    try std.testing.expectEqual(992, y20[2].im);
    try std.testing.expectEqual(-1520, y20[0].re);
    try std.testing.expectEqual(1254, y20[0].im);

    const x21 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x21);
    const y21 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y21);

    @memcpy(x21.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y21.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .RowMajor, .ConjTrans, m, n, gamma, B.ptr, n, x21.ptr, 2, delta, y21.ptr, 2);

    try std.testing.expectEqual(744, y21[0].re);
    try std.testing.expectEqual(886, y21[0].im);
    try std.testing.expectEqual(808, y21[2].re);
    try std.testing.expectEqual(972, y21[2].im);
    try std.testing.expectEqual(872, y21[4].re);
    try std.testing.expectEqual(1058, y21[4].im);
    try std.testing.expectEqual(936, y21[6].re);
    try std.testing.expectEqual(1144, y21[6].im);
    try std.testing.expectEqual(1000, y21[8].re);
    try std.testing.expectEqual(1230, y21[8].im);

    const x22 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x22);
    const y22 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y22);

    @memcpy(x22.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y22.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .RowMajor, .ConjTrans, m, n, gamma, B.ptr, n, x22.ptr, -2, delta, y22.ptr, -2);

    try std.testing.expectEqual(744, y22[8].re);
    try std.testing.expectEqual(886, y22[8].im);
    try std.testing.expectEqual(808, y22[6].re);
    try std.testing.expectEqual(972, y22[6].im);
    try std.testing.expectEqual(872, y22[4].re);
    try std.testing.expectEqual(1058, y22[4].im);
    try std.testing.expectEqual(936, y22[2].re);
    try std.testing.expectEqual(1144, y22[2].im);
    try std.testing.expectEqual(1000, y22[0].re);
    try std.testing.expectEqual(1230, y22[0].im);

    const x23 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x23);
    const y23 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y23);

    @memcpy(x23.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y23.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .ColumnMajor, .ConjTrans, m, n, gamma, B.ptr, m, x23.ptr, 2, delta, y23.ptr, 2);

    try std.testing.expectEqual(200, y23[0].re);
    try std.testing.expectEqual(246, y23[0].im);
    try std.testing.expectEqual(456, y23[2].re);
    try std.testing.expectEqual(572, y23[2].im);
    try std.testing.expectEqual(712, y23[4].re);
    try std.testing.expectEqual(898, y23[4].im);
    try std.testing.expectEqual(968, y23[6].re);
    try std.testing.expectEqual(1224, y23[6].im);
    try std.testing.expectEqual(1224, y23[8].re);
    try std.testing.expectEqual(1550, y23[8].im);

    const x24 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x24);
    const y24 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y24);

    @memcpy(x24.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y24.ptr, &[_]Complex(f64){
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

    blas.gemv(Complex(f64), .ColumnMajor, .ConjTrans, m, n, gamma, B.ptr, m, x24.ptr, -2, delta, y24.ptr, -2);

    try std.testing.expectEqual(200, y24[8].re);
    try std.testing.expectEqual(246, y24[8].im);
    try std.testing.expectEqual(456, y24[6].re);
    try std.testing.expectEqual(572, y24[6].im);
    try std.testing.expectEqual(712, y24[4].re);
    try std.testing.expectEqual(898, y24[4].im);
    try std.testing.expectEqual(968, y24[2].re);
    try std.testing.expectEqual(1224, y24[2].im);
    try std.testing.expectEqual(1224, y24[0].re);
    try std.testing.expectEqual(1550, y24[0].im);
}
