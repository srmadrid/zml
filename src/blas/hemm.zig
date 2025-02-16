const std = @import("std");
const core = @import("../core/core.zig");
const Order = @import("../ndarray/ndarray.zig").Order;
const Side = @import("../ndarray/ndarray.zig").Side;
const Uplo = @import("../ndarray/ndarray.zig").Uplo;
const blas = @import("blas.zig");

const scalar = core.supported.scalar;

pub inline fn hemm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    var SIDE = side;
    var UPLO = uplo;
    if (order == .RowMajor) {
        M = n;
        N = m;
        SIDE = if (side == .Left) .Right else .Left;
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (SIDE == .Left and lda < @max(1, M)) return;
    if (SIDE == .Right and lda < @max(1, N)) return;
    if (ldb < @max(1, M)) return;
    if (ldc < @max(1, M)) return;

    switch (supported) {
        .BuiltinBool => @compileError("blas.hemm does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("blas.hemm does not support int or float."),
        .Complex => {
            if (alpha.re == 0 and alpha.im == 0 and beta.re == 1 and beta.im == 0) return;

            if (alpha.re == 0 and alpha.im == 0) {
                if (beta.re == 0 and beta.im == 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)].re = 0;
                            C[@intCast(ijcj)].im = 0;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                } else if (beta.re != 1 or beta.im != 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            const tmp = C[@intCast(ijcj)].re * beta.re - C[@intCast(ijcj)].im * beta.im;
                            C[@intCast(ijcj)].im = C[@intCast(ijcj)].re * beta.im + C[@intCast(ijcj)].im * beta.re;
                            C[@intCast(ijcj)].re = tmp;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                }

                return;
            }

            if (SIDE == .Left) {
                if (UPLO == .Upper) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var jai: isize = 0;
                        var ibij: isize = jbj;
                        var icij: isize = jcj;
                        while (i < M) {
                            var t0: T = undefined;
                            t0.re = alpha.re * B[@intCast(ibij)].re - alpha.im * B[@intCast(ibij)].im;
                            t0.im = alpha.re * B[@intCast(ibij)].im + alpha.im * B[@intCast(ibij)].re;
                            var t1 = T.init(0, 0);

                            var k: isize = 0;
                            var iaki: isize = jai;
                            var ibkj: isize = jbj;
                            var ickj: isize = jcj;
                            while (k < i) {
                                C[@intCast(ickj)].re += t0.re * A[@intCast(iaki)].re - t0.im * A[@intCast(iaki)].im;
                                C[@intCast(ickj)].im += t0.re * A[@intCast(iaki)].im + t0.im * A[@intCast(iaki)].re;
                                t1.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                t1.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                k += 1;
                                iaki += 1;
                                ibkj += 1;
                                ickj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            k = i + jai;
                            C[@intCast(icij)].re += t0.re * A[@intCast(k)].re + t1.re * alpha.re - t1.im * alpha.im;
                            C[@intCast(icij)].im += t0.im * A[@intCast(k)].re + t1.re * alpha.im + t1.im * alpha.re;

                            i += 1;
                            jai += lda;
                            ibij += 1;
                            icij += 1;
                        }

                        j += 1;
                        jbj += ldb;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = M - 1;
                        var jai: isize = (M - 1) * lda;
                        var ibij: isize = jbj + (M - 1);
                        var icij: isize = jcj + (M - 1);
                        while (i >= 0) {
                            var t0: T = undefined;
                            t0.re = alpha.re * B[@intCast(ibij)].re - alpha.im * B[@intCast(ibij)].im;
                            t0.im = alpha.re * B[@intCast(ibij)].im + alpha.im * B[@intCast(ibij)].re;
                            var t1 = T.init(0, 0);

                            var k: isize = i + 1;
                            var iaki: isize = jai + i + 1;
                            var ibkj: isize = jbj + i + 1;
                            var ickj: isize = jcj + i + 1;
                            while (k < M) {
                                C[@intCast(ickj)].re += t0.re * A[@intCast(iaki)].re - t0.im * A[@intCast(iaki)].im;
                                C[@intCast(ickj)].im += t0.re * A[@intCast(iaki)].im + t0.im * A[@intCast(iaki)].re;
                                t1.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                t1.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                k += 1;
                                iaki += 1;
                                ibkj += 1;
                                ickj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            k = i + jai;
                            C[@intCast(icij)].re += t0.re * A[@intCast(k)].re + t1.re * alpha.re - t1.im * alpha.im;
                            C[@intCast(icij)].im += t0.im * A[@intCast(k)].re + t1.re * alpha.im + t1.im * alpha.re;

                            i -= 1;
                            jai -= lda;
                            ibij -= 1;
                            icij -= 1;
                        }

                        j += 1;
                        jbj += ldb;
                        jcj += ldc;
                    }
                }
            } else {
                if (UPLO == .Upper) {
                    var j: isize = 0;
                    var iaj: isize = 0;
                    var jaj: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = j + jaj;
                        var t0: T = undefined;
                        t0.re = alpha.re * A[@intCast(i)].re;
                        t0.im = alpha.im * A[@intCast(i)].re;

                        i = 0;
                        var ibij: isize = jbj;
                        var icij: isize = jcj;
                        while (i < M) {
                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += t0.re * B[@intCast(ibij)].re - t0.im * B[@intCast(ibij)].im;
                            C[@intCast(icij)].im += t0.re * B[@intCast(ibij)].im + t0.im * B[@intCast(ibij)].re;

                            i += 1;
                            ibij += 1;
                            icij += 1;
                        }

                        var k: isize = 0;
                        var iakj: isize = jaj;
                        var jbk: isize = 0;
                        while (k < j) {
                            t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                            t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                C[@intCast(icij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                i += 1;
                                ibik += 1;
                                icij += 1;
                            }

                            k += 1;
                            iakj += 1;
                            jbk += ldb;
                        }

                        k = j + 1;
                        var iajk: isize = iaj + (j + 1) * lda;
                        jbk = (j + 1) * ldb;
                        while (k < N) {
                            t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                            t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                C[@intCast(icij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                i += 1;
                                ibik += 1;
                                icij += 1;
                            }

                            k += 1;
                            iajk += lda;
                            jbk += ldb;
                        }

                        j += 1;
                        iaj += 1;
                        jaj += lda;
                        jbj += ldb;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var iaj: isize = 0;
                    var jaj: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = j + jaj;
                        var t0: T = undefined;
                        t0.re = alpha.re * A[@intCast(i)].re;
                        t0.im = alpha.im * A[@intCast(i)].re;

                        i = 0;
                        var ibij: isize = jbj;
                        var icij: isize = jcj;
                        while (i < M) {
                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += t0.re * B[@intCast(ibij)].re - t0.im * B[@intCast(ibij)].im;
                            C[@intCast(icij)].im += t0.re * B[@intCast(ibij)].im + t0.im * B[@intCast(ibij)].re;

                            i += 1;
                            ibij += 1;
                            icij += 1;
                        }

                        var k: isize = 0;
                        var iajk: isize = iaj;
                        var jbk: isize = 0;
                        while (k < j) {
                            t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                            t0.im = -alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                C[@intCast(icij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                i += 1;
                                ibik += 1;
                                icij += 1;
                            }

                            k += 1;
                            iajk += lda;
                            jbk += ldb;
                        }

                        k = j + 1;
                        var iakj: isize = j + 1 + jaj;
                        jbk = (j + 1) * ldb;
                        while (k < N) {
                            t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                            t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                C[@intCast(icij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                i += 1;
                                ibik += 1;
                                icij += 1;
                            }

                            k += 1;
                            iakj += 1;
                            jbk += ldb;
                        }

                        j += 1;
                        iaj += 1;
                        jaj += lda;
                        jbj += ldb;
                        jcj += ldc;
                    }
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.hemm only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "hemm" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const alpha = Complex(f64).init(2, 2);
    const beta = Complex(f64).init(3, 3);

    const A = try a.alloc(Complex(f64), m * m);
    defer a.free(A);
    const B = try a.alloc(Complex(f64), n * n);
    defer a.free(B);
    const C = try a.alloc(Complex(f64), m * n);
    defer a.free(C);
    const D = try a.alloc(Complex(f64), m * n);
    defer a.free(D);

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
    });
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
        Complex(f64).init(21, 21),
        Complex(f64).init(22, 22),
        Complex(f64).init(23, 23),
        Complex(f64).init(24, 24),
        Complex(f64).init(25, 25),
    });
    @memcpy(C.ptr, &[_]Complex(f64){
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
    @memcpy(D.ptr, &[_]Complex(f64){
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

    blas.hemm(Complex(f64), .RowMajor, .Left, .Upper, m, n, alpha, A.ptr, m, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(-436, D[0].re);
    try std.testing.expectEqual(446, D[0].im);
    try std.testing.expectEqual(-472, D[1].re);
    try std.testing.expectEqual(492, D[1].im);
    try std.testing.expectEqual(-508, D[2].re);
    try std.testing.expectEqual(538, D[2].im);
    try std.testing.expectEqual(-544, D[3].re);
    try std.testing.expectEqual(584, D[3].im);
    try std.testing.expectEqual(-580, D[4].re);
    try std.testing.expectEqual(630, D[4].im);
    try std.testing.expectEqual(-812, D[5].re);
    try std.testing.expectEqual(1008, D[5].im);
    try std.testing.expectEqual(-864, D[6].re);
    try std.testing.expectEqual(1106, D[6].im);
    try std.testing.expectEqual(-916, D[7].re);
    try std.testing.expectEqual(1204, D[7].im);
    try std.testing.expectEqual(-968, D[8].re);
    try std.testing.expectEqual(1302, D[8].im);
    try std.testing.expectEqual(-1020, D[9].re);
    try std.testing.expectEqual(1400, D[9].im);
    try std.testing.expectEqual(-588, D[10].re);
    try std.testing.expectEqual(1498, D[10].im);
    try std.testing.expectEqual(-596, D[11].re);
    try std.testing.expectEqual(1636, D[11].im);
    try std.testing.expectEqual(-604, D[12].re);
    try std.testing.expectEqual(1774, D[12].im);
    try std.testing.expectEqual(-612, D[13].re);
    try std.testing.expectEqual(1912, D[13].im);
    try std.testing.expectEqual(-620, D[14].re);
    try std.testing.expectEqual(2050, D[14].im);
    try std.testing.expectEqual(736, D[15].re);
    try std.testing.expectEqual(1856, D[15].im);
    try std.testing.expectEqual(832, D[16].re);
    try std.testing.expectEqual(2022, D[16].im);
    try std.testing.expectEqual(928, D[17].re);
    try std.testing.expectEqual(2188, D[17].im);
    try std.testing.expectEqual(1024, D[18].re);
    try std.testing.expectEqual(2354, D[18].im);
    try std.testing.expectEqual(1120, D[19].re);
    try std.testing.expectEqual(2520, D[19].im);

    blas.hemm(Complex(f64), .ColumnMajor, .Left, .Upper, m, n, alpha, A.ptr, m, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(-3002, D[0].re);
    try std.testing.expectEqual(390, D[0].im);
    try std.testing.expectEqual(-3216, D[1].re);
    try std.testing.expectEqual(472, D[1].im);
    try std.testing.expectEqual(-3262, D[2].re);
    try std.testing.expectEqual(578, D[2].im);
    try std.testing.expectEqual(-3040, D[3].re);
    try std.testing.expectEqual(720, D[3].im);
    try std.testing.expectEqual(-4418, D[4].re);
    try std.testing.expectEqual(958, D[4].im);
    try std.testing.expectEqual(-6088, D[5].re);
    try std.testing.expectEqual(1560, D[5].im);
    try std.testing.expectEqual(-5970, D[6].re);
    try std.testing.expectEqual(1934, D[6].im);
    try std.testing.expectEqual(-5344, D[7].re);
    try std.testing.expectEqual(2392, D[7].im);
    try std.testing.expectEqual(-8030, D[8].re);
    try std.testing.expectEqual(2258, D[8].im);
    try std.testing.expectEqual(-8192, D[9].re);
    try std.testing.expectEqual(2672, D[9].im);
    try std.testing.expectEqual(-6254, D[10].re);
    try std.testing.expectEqual(4658, D[10].im);
    try std.testing.expectEqual(-5008, D[11].re);
    try std.testing.expectEqual(5576, D[11].im);
    try std.testing.expectEqual(-8786, D[12].re);
    try std.testing.expectEqual(5214, D[12].im);
    try std.testing.expectEqual(-8808, D[13].re);
    try std.testing.expectEqual(5992, D[13].im);
    try std.testing.expectEqual(-7942, D[14].re);
    try std.testing.expectEqual(6938, D[14].im);
    try std.testing.expectEqual(-1000, D[15].re);
    try std.testing.expectEqual(11160, D[15].im);
    try std.testing.expectEqual(-5654, D[16].re);
    try std.testing.expectEqual(10714, D[16].im);
    try std.testing.expectEqual(-5320, D[17].re);
    try std.testing.expectEqual(12000, D[17].im);
    try std.testing.expectEqual(-3858, D[18].re);
    try std.testing.expectEqual(13502, D[18].im);
    try std.testing.expectEqual(-1168, D[19].re);
    try std.testing.expectEqual(15232, D[19].im);

    blas.hemm(Complex(f64), .RowMajor, .Left, .Lower, m, n, alpha, A.ptr, m, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(-8828, D[0].re);
    try std.testing.expectEqual(-6484, D[0].im);
    try std.testing.expectEqual(-9608, D[1].re);
    try std.testing.expectEqual(-6768, D[1].im);
    try std.testing.expectEqual(-9956, D[2].re);
    try std.testing.expectEqual(-6476, D[2].im);
    try std.testing.expectEqual(-9608, D[3].re);
    try std.testing.expectEqual(-5272, D[3].im);
    try std.testing.expectEqual(-14348, D[4].re);
    try std.testing.expectEqual(-8580, D[4].im);
    try std.testing.expectEqual(-21628, D[5].re);
    try std.testing.expectEqual(-12084, D[5].im);
    try std.testing.expectEqual(-22320, D[6].re);
    try std.testing.expectEqual(-10468, D[6].im);
    try std.testing.expectEqual(-21740, D[7].re);
    try std.testing.expectEqual(-7076, D[7].im);
    try std.testing.expectEqual(-29320, D[8].re);
    try std.testing.expectEqual(-15396, D[8].im);
    try std.testing.expectEqual(-30972, D[9].re);
    try std.testing.expectEqual(-14500, D[9].im);
    try std.testing.expectEqual(-32052, D[10].re);
    try std.testing.expectEqual(-3068, D[10].im);
    try std.testing.expectEqual(-31084, D[11].re);
    try std.testing.expectEqual(3604, D[11].im);
    try std.testing.expectEqual(-41348, D[12].re);
    try std.testing.expectEqual(-8636, D[12].im);
    try std.testing.expectEqual(-43764, D[13].re);
    try std.testing.expectEqual(-6188, D[13].im);
    try std.testing.expectEqual(-44020, D[14].re);
    try std.testing.expectEqual(-572, D[14].im);
    try std.testing.expectEqual(-37528, D[15].re);
    try std.testing.expectEqual(32552, D[15].im);
    try std.testing.expectEqual(-50320, D[16].re);
    try std.testing.expectEqual(17484, D[16].im);
    try std.testing.expectEqual(-53344, D[17].re);
    try std.testing.expectEqual(22576, D[17].im);
    try std.testing.expectEqual(-53632, D[18].re);
    try std.testing.expectEqual(31700, D[18].im);
    try std.testing.expectEqual(-50920, D[19].re);
    try std.testing.expectEqual(45192, D[19].im);

    blas.hemm(Complex(f64), .ColumnMajor, .Left, .Lower, m, n, alpha, A.ptr, m, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(-6916, D[0].re);
    try std.testing.expectEqual(-45816, D[0].im);
    try std.testing.expectEqual(-8316, D[1].re);
    try std.testing.expectEqual(-48860, D[1].im);
    try std.testing.expectEqual(-10316, D[2].re);
    try std.testing.expectEqual(-48904, D[2].im);
    try std.testing.expectEqual(-13232, D[3].re);
    try std.testing.expectEqual(-44160, D[3].im);
    try std.testing.expectEqual(-17044, D[4].re);
    try std.testing.expectEqual(-68504, D[4].im);
    try std.testing.expectEqual(-28220, D[5].re);
    try std.testing.expectEqual(-100500, D[5].im);
    try std.testing.expectEqual(-35400, D[6].re);
    try std.testing.expectEqual(-97444, D[6].im);
    try std.testing.expectEqual(-44600, D[7].re);
    try std.testing.expectEqual(-85328, D[7].im);
    try std.testing.expectEqual(-41368, D[8].re);
    try std.testing.expectEqual(-133708, D[8].im);
    try std.testing.expectEqual(-48796, D[9].re);
    try std.testing.expectEqual(-135412, D[9].im);
    try std.testing.expectEqual(-86764, D[10].re);
    try std.testing.expectEqual(-103912, D[10].im);
    try std.testing.expectEqual(-105056, D[11].re);
    try std.testing.expectEqual(-80680, D[11].im);
    try std.testing.expectEqual(-97588, D[12].re);
    try std.testing.expectEqual(-149352, D[12].im);
    try std.testing.expectEqual(-111900, D[13].re);
    try std.testing.expectEqual(-148484, D[13].im);
    try std.testing.expectEqual(-130124, D[14].re);
    try std.testing.expectEqual(-131800, D[14].im);
    try std.testing.expectEqual(-211616, D[15].re);
    try std.testing.expectEqual(-12528, D[15].im);
    try std.testing.expectEqual(-202720, D[16].re);
    try std.testing.expectEqual(-97748, D[16].im);
    try std.testing.expectEqual(-226724, D[17].re);
    try std.testing.expectEqual(-90564, D[17].im);
    try std.testing.expectEqual(-255744, D[18].re);
    try std.testing.expectEqual(-63292, D[18].im);
    try std.testing.expectEqual(-290096, D[19].re);
    try std.testing.expectEqual(-14144, D[19].im);

    blas.hemm(Complex(f64), .RowMajor, .Right, .Upper, m, n, alpha, B.ptr, n, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(116916, D[0].re);
    try std.testing.expectEqual(-157976, D[0].im);
    try std.testing.expectEqual(122064, D[1].re);
    try std.testing.expectEqual(-171024, D[1].im);
    try std.testing.expectEqual(116212, D[2].re);
    try std.testing.expectEqual(-176904, D[2].im);
    try std.testing.expectEqual(92928, D[3].re);
    try std.testing.expectEqual(-171216, D[3].im);
    try std.testing.expectEqual(153780, D[4].re);
    try std.testing.expectEqual(-255544, D[4].im);
    try std.testing.expectEqual(217336, D[5].re);
    try std.testing.expectEqual(-385640, D[5].im);
    try std.testing.expectEqual(187064, D[6].re);
    try std.testing.expectEqual(-397308, D[6].im);
    try std.testing.expectEqual(122992, D[7].re);
    try std.testing.expectEqual(-387968, D[7].im);
    try std.testing.expectEqual(277024, D[8].re);
    try std.testing.expectEqual(-522948, D[8].im);
    try std.testing.expectEqual(258248, D[9].re);
    try std.testing.expectEqual(-550024, D[9].im);
    try std.testing.expectEqual(52220, D[10].re);
    try std.testing.expectEqual(-571208, D[10].im);
    try std.testing.expectEqual(-71696, D[11].re);
    try std.testing.expectEqual(-555264, D[11].im);
    try std.testing.expectEqual(156460, D[12].re);
    try std.testing.expectEqual(-737944, D[12].im);
    try std.testing.expectEqual(109616, D[13].re);
    try std.testing.expectEqual(-777552, D[13].im);
    try std.testing.expectEqual(2428, D[14].re);
    try std.testing.expectEqual(-781672, D[14].im);
    try std.testing.expectEqual(-596208, D[15].re);
    try std.testing.expectEqual(-671312, D[15].im);
    try std.testing.expectEqual(-312984, D[16].re);
    try std.testing.expectEqual(-898740, D[16].im);
    try std.testing.expectEqual(-406952, D[17].re);
    try std.testing.expectEqual(-947928, D[17].im);
    try std.testing.expectEqual(-577632, D[18].re);
    try std.testing.expectEqual(-952188, D[18].im);
    try std.testing.expectEqual(-831456, D[19].re);
    try std.testing.expectEqual(-907120, D[19].im);

    blas.hemm(Complex(f64), .ColumnMajor, .Right, .Upper, m, n, alpha, B.ptr, n, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(827452, D[0].re);
    try std.testing.expectEqual(-120400, D[0].im);
    try std.testing.expectEqual(882256, D[1].re);
    try std.testing.expectEqual(-143880, D[1].im);
    try std.testing.expectEqual(882556, D[2].re);
    try std.testing.expectEqual(-178856, D[2].im);
    try std.testing.expectEqual(795856, D[3].re);
    try std.testing.expectEqual(-231424, D[3].im);
    try std.testing.expectEqual(1230760, D[4].re);
    try std.testing.expectEqual(-302316, D[4].im);
    try std.testing.expectEqual(1811896, D[5].re);
    try std.testing.expectEqual(-501680, D[5].im);
    try std.testing.expectEqual(1756264, D[6].re);
    try std.testing.expectEqual(-627244, D[6].im);
    try std.testing.expectEqual(1536208, D[7].re);
    try std.testing.expectEqual(-791184, D[7].im);
    try std.testing.expectEqual(2402132, D[8].re);
    try std.testing.expectEqual(-734520, D[8].im);
    try std.testing.expectEqual(2427104, D[9].re);
    try std.testing.expectEqual(-871768, D[9].im);
    try std.testing.expectEqual(1872644, D[10].re);
    try std.testing.expectEqual(-1553096, D[10].im);
    try std.testing.expectEqual(1453136, D[11].re);
    try std.testing.expectEqual(-1876704, D[11].im);
    try std.testing.expectEqual(2683792, D[12].re);
    try std.testing.expectEqual(-1740780, D[12].im);
    try std.testing.expectEqual(2661976, D[13].re);
    try std.testing.expectEqual(-1999760, D[13].im);
    try std.testing.expectEqual(2352664, D[14].re);
    try std.testing.expectEqual(-2333308, D[14].im);
    try std.testing.expectEqual(225568, D[15].re);
    try std.testing.expectEqual(-3797760, D[15].im);
    try std.testing.expectEqual(1754668, D[16].re);
    try std.testing.expectEqual(-3630872, D[16].im);
    try std.testing.expectEqual(1619968, D[17].re);
    try std.testing.expectEqual(-4059880, D[17].im);
    try std.testing.expectEqual(1120348, D[18].re);
    try std.testing.expectEqual(-4584240, D[18].im);
    try std.testing.expectEqual(223312, D[19].re);
    try std.testing.expectEqual(-5210048, D[19].im);

    blas.hemm(Complex(f64), .RowMajor, .Right, .Lower, m, n, alpha, B.ptr, n, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(2842700, D[0].re);
    try std.testing.expectEqual(2122016, D[0].im);
    try std.testing.expectEqual(3077576, D[1].re);
    try std.testing.expectEqual(2216064, D[1].im);
    try std.testing.expectEqual(3183628, D[2].re);
    try std.testing.expectEqual(2112144, D[2].im);
    try std.testing.expectEqual(3081776, D[3].re);
    try std.testing.expectEqual(1694496, D[3].im);
    try std.testing.expectEqual(4600148, D[4].re);
    try std.testing.expectEqual(2786752, D[4].im);
    try std.testing.expectEqual(6938792, D[5].re);
    try std.testing.expectEqual(3932608, D[5].im);
    try std.testing.expectEqual(7148792, D[6].re);
    try std.testing.expectEqual(3389276, D[6].im);
    try std.testing.expectEqual(6981208, D[7].re);
    try std.testing.expectEqual(2237656, D[7].im);
    try std.testing.expectEqual(9410432, D[8].re);
    try std.testing.expectEqual(5005916, D[8].im);
    try std.testing.expectEqual(9899336, D[9].re);
    try std.testing.expectEqual(4669728, D[9].im);
    try std.testing.expectEqual(10274204, D[10].re);
    try std.testing.expectEqual(961704, D[10].im);
    try std.testing.expectEqual(9986888, D[11].re);
    try std.testing.expectEqual(-1267208, D[11].im);
    try std.testing.expectEqual(13272388, D[12].re);
    try std.testing.expectEqual(2833160, D[12].im);
    try std.testing.expectEqual(13986224, D[13].re);
    try std.testing.expectEqual(1991608, D[13].im);
    try std.testing.expectEqual(14062436, D[14].re);
    try std.testing.expectEqual(64088, D[14].im);
    try std.testing.expectEqual(12065888, D[15].re);
    try std.testing.expectEqual(-10712416, D[15].im);
    try std.testing.expectEqual(16153088, D[16].re);
    try std.testing.expectEqual(-5623836, D[16].im);
    try std.testing.expectEqual(17037856, D[17].re);
    try std.testing.expectEqual(-7314072, D[17].im);
    try std.testing.expectEqual(17115320, D[18].re);
    try std.testing.expectEqual(-10384836, D[18].im);
    try std.testing.expectEqual(16306400, D[19].re);
    try std.testing.expectEqual(-14951888, D[19].im);

    blas.hemm(Complex(f64), .ColumnMajor, .Right, .Lower, m, n, alpha, B.ptr, n, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(2161356, D[0].re);
    try std.testing.expectEqual(14894848, D[0].im);
    try std.testing.expectEqual(2583784, D[1].re);
    try std.testing.expectEqual(15881680, D[1].im);
    try std.testing.expectEqual(3213644, D[2].re);
    try std.testing.expectEqual(15888136, D[2].im);
    try std.testing.expectEqual(4160976, D[3].re);
    try std.testing.expectEqual(14329696, D[3].im);
    try std.testing.expectEqual(5438760, D[4].re);
    try std.testing.expectEqual(22162284, D[4].im);
    try std.testing.expectEqual(9017024, D[5].re);
    try std.testing.expectEqual(32615928, D[5].im);
    try std.testing.expectEqual(11276920, D[6].re);
    try std.testing.expectEqual(31616076, D[6].im);
    try std.testing.expectEqual(14228928, D[7].re);
    try std.testing.expectEqual(27658608, D[7].im);
    try std.testing.expectEqual(13211972, D[8].re);
    try std.testing.expectEqual(43251432, D[8].im);
    try std.testing.expectEqual(15687176, D[9].re);
    try std.testing.expectEqual(43709792, D[9].im);
    try std.testing.expectEqual(27935780, D[10].re);
    try std.testing.expectEqual(33710536, D[10].im);
    try std.testing.expectEqual(33760496, D[11].re);
    try std.testing.expectEqual(26162064, D[11].im);
    try std.testing.expectEqual(31317024, D[12].re);
    try std.testing.expectEqual(48319692, D[12].im);
    try std.testing.expectEqual(35983216, D[13].re);
    try std.testing.expectEqual(47936808, D[13].im);
    try std.testing.expectEqual(41994440, D[14].re);
    try std.testing.expectEqual(42383148, D[14].im);
    try std.testing.expectEqual(68334336, D[15].re);
    try std.testing.expectEqual(4064256, D[15].im);
    try std.testing.expectEqual(65332572, D[16].re);
    try std.testing.expectEqual(31591256, D[16].im);
    try std.testing.expectEqual(73057784, D[17].re);
    try std.testing.expectEqual(29175152, D[17].im);
    try std.testing.expectEqual(82502668, D[18].re);
    try std.testing.expectEqual(20195552, D[18].im);
    try std.testing.expectEqual(93777264, D[19].re);
    try std.testing.expectEqual(4067936, D[19].im);
}
