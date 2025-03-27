const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

const Numeric = core.types.Numeric;

pub inline fn herk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Numeric(T), A: [*]const T, lda: isize, beta: Numeric(T), C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or ((alpha == 0 or k <= 0) and beta == 1) or trans == .Trans or trans == .ConjNoTrans) return;

    var UPLO = uplo;
    var TRANS = trans;
    var NROWA = if (trans == .NoTrans) n else k;
    const ldcp1 = ldc + 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANS = if (trans == .NoTrans) .ConjTrans else .NoTrans;
        NROWA = if (trans == .NoTrans) k else n;
    }

    if (lda < @max(1, NROWA)) return;
    if (ldc < @max(1, n)) return;

    switch (numericType) {
        .bool => @compileError("blas.herk does not support bool."),
        .int, .float => @compileError("blas.herk does not support int or float."),
        .cfloat => {
            if (alpha == 0) {
                if (UPLO == .Upper) {
                    if (beta == 0) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = 0;
                            var icij: isize = jcj;
                            while (i <= j) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;

                                i += 1;
                                icij += 1;
                            }

                            j += 1;
                            jcj += ldc;
                        }
                    } else if (beta != 1) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = 0;
                            var icij: isize = jcj;
                            while (i < j) {
                                C[@intCast(icij)].re *= beta;
                                C[@intCast(icij)].im *= beta;

                                i += 1;
                                icij += 1;
                            }

                            C[@intCast(icij)].re *= beta;
                            C[@intCast(icij)].im = 0;

                            j += 1;
                            jcj += ldc;
                        }
                    }
                } else {
                    if (beta == 0) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = j;
                            var icij: isize = jcj;
                            while (i < n) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;

                                i += 1;
                                icij += 1;
                            }

                            j += 1;
                            jcj += ldcp1;
                        }
                    } else if (beta != 1) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            C[@intCast(jcj)].re *= beta;
                            C[@intCast(jcj)].im = 0;

                            var i: isize = j + 1;
                            var icij: isize = jcj + 1;
                            while (i < n) {
                                C[@intCast(icij)].re *= beta;
                                C[@intCast(icij)].im *= beta;

                                i += 1;
                                icij += 1;
                            }

                            j += 1;
                            jcj += ldcp1;
                        }
                    }
                }

                return;
            }

            if (UPLO == .Upper) {
                if (TRANS == .NoTrans) {
                    var j: isize = 0;
                    var iaj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        if (beta == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < j) {
                                Cpjcj[@intCast(icij)].re = 0;
                                Cpjcj[@intCast(icij)].im = 0;

                                icij += 1;
                            }
                        } else if (beta != 1) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < j) {
                                Cpjcj[@intCast(icij)].re *= beta;
                                Cpjcj[@intCast(icij)].im *= beta;

                                icij += 1;
                            }
                        }

                        var icij: isize = j + jcj;

                        if (beta == 0) {
                            C[@intCast(icij)].re = 0;
                        } else if (beta != 1) {
                            C[@intCast(icij)].re *= beta;
                        }
                        C[@intCast(icij)].im = 0;

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var jal: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha * A[@intCast(iajl)].re;
                            t0.im = -alpha * A[@intCast(iajl)].im;

                            var i: isize = 0;
                            var iail: isize = jal;
                            icij = jcj;
                            while (i < j) {
                                C[@intCast(icij)].re += t0.re * A[@intCast(iail)].re - t0.im * A[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.re * A[@intCast(iail)].im + t0.im * A[@intCast(iail)].re;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            C[@intCast(icij)].re += t0.re * A[@intCast(iajl)].re - t0.im * A[@intCast(iajl)].im;
                            C[@intCast(icij)].im = 0;

                            l += 1;
                            iajl += lda;
                            jal += lda;
                        }

                        j += 1;
                        iaj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jaj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        var i: isize = 0;
                        var jai: isize = 0;
                        var icij: isize = jcj;
                        while (i < j) {
                            var t0 = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re + A[@intCast(iali)].im * A[@intCast(ialj)].im;
                                t0.im += A[@intCast(iali)].re * A[@intCast(ialj)].im - A[@intCast(iali)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)].re *= beta;
                                C[@intCast(icij)].im *= beta;
                            }

                            C[@intCast(icij)].re += alpha * t0.re;
                            C[@intCast(icij)].im += alpha * t0.im;

                            i += 1;
                            jai += lda;
                            icij += 1;
                        }

                        var t0 = T.init(0, 0);

                        var l: isize = 0;
                        var iali: isize = jai;
                        var ialj: isize = jaj;
                        while (l < k) {
                            t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re + A[@intCast(iali)].im * A[@intCast(ialj)].im;

                            l += 1;
                            iali += 1;
                            ialj += 1;
                        }

                        if (beta == 0) {
                            C[@intCast(icij)].re = 0;
                        } else if (beta != 1) {
                            C[@intCast(icij)].re *= beta;
                        }

                        C[@intCast(icij)].re += alpha * t0.re;
                        C[@intCast(icij)].im = 0;

                        j += 1;
                        jaj += lda;
                        jcj += ldc;
                    }
                }
            } else {
                if (TRANS == .NoTrans) {
                    var j: isize = 0;
                    var iaj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        var icij: isize = j + jcj;

                        if (beta == 0) {
                            C[@intCast(icij)].re = 0;
                            C[@intCast(icij)].im = 0;

                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(icij + 1)]);
                            while (icj < n - j - 1) {
                                Cpicij[@intCast(icj)].re = 0;
                                Cpicij[@intCast(icj)].im = 0;

                                icj += 1;
                            }
                        } else if (beta != 1) {
                            C[@intCast(icij)].re *= beta;
                            C[@intCast(icij)].im *= beta;

                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(icij + 1)]);
                            while (icj < n - j - 1) {
                                Cpicij[@intCast(icj)].re *= beta;
                                Cpicij[@intCast(icj)].im *= beta;

                                icj += 1;
                            }
                        }

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var jal: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha * A[@intCast(iajl)].re;
                            t0.im = -alpha * A[@intCast(iajl)].im;

                            var iail: isize = j + jal;
                            icij = j + jcj;

                            C[@intCast(icij)].re += t0.re * A[@intCast(iajl)].re - t0.im * A[@intCast(iajl)].im;
                            C[@intCast(icij)].im = 0;

                            iail += 1;
                            icij += 1;

                            var i: isize = j + 1;
                            while (i < n) {
                                C[@intCast(icij)].re += t0.re * A[@intCast(iail)].re - t0.im * A[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.re * A[@intCast(iail)].im + t0.im * A[@intCast(iail)].re;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            iajl += lda;
                            jal += lda;
                        }

                        j += 1;
                        iaj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jaj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        var jai: isize = j * lda;
                        var icij: isize = j + jcj;

                        var t0 = T.init(0, 0);

                        var l: isize = 0;
                        var iali: isize = jai;
                        var ialj: isize = jaj;
                        while (l < k) {
                            t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re + A[@intCast(iali)].im * A[@intCast(ialj)].im;

                            l += 1;
                            iali += 1;
                            ialj += 1;
                        }

                        if (beta == 0) {
                            C[@intCast(icij)].re = 0;
                        } else if (beta != 1) {
                            C[@intCast(icij)].re *= beta;
                        }

                        C[@intCast(icij)].re += alpha * t0.re;
                        C[@intCast(icij)].im = 0;

                        jai += lda;
                        icij += 1;

                        var i: isize = j + 1;
                        while (i < n) {
                            t0 = T.init(0, 0);

                            l = 0;
                            iali = jai;
                            ialj = jaj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re + A[@intCast(iali)].im * A[@intCast(ialj)].im;
                                t0.im += A[@intCast(iali)].re * A[@intCast(ialj)].im - A[@intCast(iali)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)].re *= beta;
                                C[@intCast(icij)].im *= beta;
                            }

                            C[@intCast(icij)].re += alpha * t0.re;
                            C[@intCast(icij)].im += alpha * t0.im;

                            i += 1;
                            jai += lda;
                            icij += 1;
                        }

                        j += 1;
                        jaj += lda;
                        jcj += ldc;
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.herk only supports simple types."),
        .unsupported => unreachable,
    }
}

test herk {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 3;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(Complex(f64), n * k);
    defer a.free(A);
    const B = try a.alloc(Complex(f64), n * n);
    defer a.free(B);

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

    blas.herk(Complex(f64), .RowMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(59, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(134, B[1].re);
    try std.testing.expectEqual(6, B[1].im);
    try std.testing.expectEqual(209, B[2].re);
    try std.testing.expectEqual(9, B[2].im);
    try std.testing.expectEqual(284, B[3].re);
    try std.testing.expectEqual(12, B[3].im);
    try std.testing.expectEqual(359, B[4].re);
    try std.testing.expectEqual(15, B[4].im);
    try std.testing.expectEqual(6, B[5].re);
    try std.testing.expectEqual(6, B[5].im);
    try std.testing.expectEqual(329, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(512, B[7].re);
    try std.testing.expectEqual(24, B[7].im);
    try std.testing.expectEqual(695, B[8].re);
    try std.testing.expectEqual(27, B[8].im);
    try std.testing.expectEqual(878, B[9].re);
    try std.testing.expectEqual(30, B[9].im);
    try std.testing.expectEqual(11, B[10].re);
    try std.testing.expectEqual(11, B[10].im);
    try std.testing.expectEqual(12, B[11].re);
    try std.testing.expectEqual(12, B[11].im);
    try std.testing.expectEqual(815, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(1106, B[13].re);
    try std.testing.expectEqual(42, B[13].im);
    try std.testing.expectEqual(1397, B[14].re);
    try std.testing.expectEqual(45, B[14].im);
    try std.testing.expectEqual(16, B[15].re);
    try std.testing.expectEqual(16, B[15].im);
    try std.testing.expectEqual(17, B[16].re);
    try std.testing.expectEqual(17, B[16].im);
    try std.testing.expectEqual(18, B[17].re);
    try std.testing.expectEqual(18, B[17].im);
    try std.testing.expectEqual(1517, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(1916, B[19].re);
    try std.testing.expectEqual(60, B[19].im);
    try std.testing.expectEqual(21, B[20].re);
    try std.testing.expectEqual(21, B[20].im);
    try std.testing.expectEqual(22, B[21].re);
    try std.testing.expectEqual(22, B[21].im);
    try std.testing.expectEqual(23, B[22].re);
    try std.testing.expectEqual(23, B[22].im);
    try std.testing.expectEqual(24, B[23].re);
    try std.testing.expectEqual(24, B[23].im);
    try std.testing.expectEqual(2435, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .ColumnMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(809, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(134, B[1].re);
    try std.testing.expectEqual(6, B[1].im);
    try std.testing.expectEqual(209, B[2].re);
    try std.testing.expectEqual(9, B[2].im);
    try std.testing.expectEqual(284, B[3].re);
    try std.testing.expectEqual(12, B[3].im);
    try std.testing.expectEqual(359, B[4].re);
    try std.testing.expectEqual(15, B[4].im);
    try std.testing.expectEqual(722, B[5].re);
    try std.testing.expectEqual(18, B[5].im);
    try std.testing.expectEqual(1775, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(512, B[7].re);
    try std.testing.expectEqual(24, B[7].im);
    try std.testing.expectEqual(695, B[8].re);
    try std.testing.expectEqual(27, B[8].im);
    try std.testing.expectEqual(878, B[9].re);
    try std.testing.expectEqual(30, B[9].im);
    try std.testing.expectEqual(809, B[10].re);
    try std.testing.expectEqual(33, B[10].im);
    try std.testing.expectEqual(908, B[11].re);
    try std.testing.expectEqual(36, B[11].im);
    try std.testing.expectEqual(3413, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(1106, B[13].re);
    try std.testing.expectEqual(42, B[13].im);
    try std.testing.expectEqual(1397, B[14].re);
    try std.testing.expectEqual(45, B[14].im);
    try std.testing.expectEqual(896, B[15].re);
    try std.testing.expectEqual(48, B[15].im);
    try std.testing.expectEqual(1007, B[16].re);
    try std.testing.expectEqual(51, B[16].im);
    try std.testing.expectEqual(1118, B[17].re);
    try std.testing.expectEqual(54, B[17].im);
    try std.testing.expectEqual(5723, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(1916, B[19].re);
    try std.testing.expectEqual(60, B[19].im);
    try std.testing.expectEqual(983, B[20].re);
    try std.testing.expectEqual(63, B[20].im);
    try std.testing.expectEqual(1106, B[21].re);
    try std.testing.expectEqual(66, B[21].im);
    try std.testing.expectEqual(1229, B[22].re);
    try std.testing.expectEqual(69, B[22].im);
    try std.testing.expectEqual(1352, B[23].re);
    try std.testing.expectEqual(72, B[23].im);
    try std.testing.expectEqual(8705, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .RowMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(3059, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(722, B[5].re);
    try std.testing.expectEqual(18, B[5].im);
    try std.testing.expectEqual(6113, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(809, B[10].re);
    try std.testing.expectEqual(33, B[10].im);
    try std.testing.expectEqual(908, B[11].re);
    try std.testing.expectEqual(36, B[11].im);
    try std.testing.expectEqual(11207, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(896, B[15].re);
    try std.testing.expectEqual(48, B[15].im);
    try std.testing.expectEqual(1007, B[16].re);
    try std.testing.expectEqual(51, B[16].im);
    try std.testing.expectEqual(1118, B[17].re);
    try std.testing.expectEqual(54, B[17].im);
    try std.testing.expectEqual(18341, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(983, B[20].re);
    try std.testing.expectEqual(63, B[20].im);
    try std.testing.expectEqual(1106, B[21].re);
    try std.testing.expectEqual(66, B[21].im);
    try std.testing.expectEqual(1229, B[22].re);
    try std.testing.expectEqual(69, B[22].im);
    try std.testing.expectEqual(1352, B[23].re);
    try std.testing.expectEqual(72, B[23].im);
    try std.testing.expectEqual(27515, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(9233, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(2294, B[5].re);
    try std.testing.expectEqual(54, B[5].im);
    try std.testing.expectEqual(18647, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(2627, B[10].re);
    try std.testing.expectEqual(99, B[10].im);
    try std.testing.expectEqual(3212, B[11].re);
    try std.testing.expectEqual(108, B[11].im);
    try std.testing.expectEqual(34397, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(2960, B[15].re);
    try std.testing.expectEqual(144, B[15].im);
    try std.testing.expectEqual(3689, B[16].re);
    try std.testing.expectEqual(153, B[16].im);
    try std.testing.expectEqual(4418, B[17].re);
    try std.testing.expectEqual(162, B[17].im);
    try std.testing.expectEqual(56483, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(3293, B[20].re);
    try std.testing.expectEqual(189, B[20].im);
    try std.testing.expectEqual(4166, B[21].re);
    try std.testing.expectEqual(198, B[21].im);
    try std.testing.expectEqual(5039, B[22].re);
    try std.testing.expectEqual(207, B[22].im);
    try std.testing.expectEqual(5912, B[23].re);
    try std.testing.expectEqual(216, B[23].im);
    try std.testing.expectEqual(84905, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .RowMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(27755, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(7010, B[5].re);
    try std.testing.expectEqual(162, B[5].im);
    try std.testing.expectEqual(56249, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(8081, B[10].re);
    try std.testing.expectEqual(297, B[10].im);
    try std.testing.expectEqual(10124, B[11].re);
    try std.testing.expectEqual(324, B[11].im);
    try std.testing.expectEqual(103967, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(9152, B[15].re);
    try std.testing.expectEqual(432, B[15].im);
    try std.testing.expectEqual(11735, B[16].re);
    try std.testing.expectEqual(459, B[16].im);
    try std.testing.expectEqual(14318, B[17].re);
    try std.testing.expectEqual(486, B[17].im);
    try std.testing.expectEqual(170909, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(10223, B[20].re);
    try std.testing.expectEqual(567, B[20].im);
    try std.testing.expectEqual(13346, B[21].re);
    try std.testing.expectEqual(594, B[21].im);
    try std.testing.expectEqual(16469, B[22].re);
    try std.testing.expectEqual(621, B[22].im);
    try std.testing.expectEqual(19592, B[23].re);
    try std.testing.expectEqual(648, B[23].im);
    try std.testing.expectEqual(257075, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .ColumnMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(83897, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(4022, B[1].re);
    try std.testing.expectEqual(54, B[1].im);
    try std.testing.expectEqual(4985, B[2].re);
    try std.testing.expectEqual(81, B[2].im);
    try std.testing.expectEqual(5948, B[3].re);
    try std.testing.expectEqual(108, B[3].im);
    try std.testing.expectEqual(6911, B[4].re);
    try std.testing.expectEqual(135, B[4].im);
    try std.testing.expectEqual(7010, B[5].re);
    try std.testing.expectEqual(162, B[5].im);
    try std.testing.expectEqual(169535, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(8096, B[7].re);
    try std.testing.expectEqual(216, B[7].im);
    try std.testing.expectEqual(10079, B[8].re);
    try std.testing.expectEqual(243, B[8].im);
    try std.testing.expectEqual(12062, B[9].re);
    try std.testing.expectEqual(270, B[9].im);
    try std.testing.expectEqual(8081, B[10].re);
    try std.testing.expectEqual(297, B[10].im);
    try std.testing.expectEqual(10124, B[11].re);
    try std.testing.expectEqual(324, B[11].im);
    try std.testing.expectEqual(312869, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(14210, B[13].re);
    try std.testing.expectEqual(378, B[13].im);
    try std.testing.expectEqual(17213, B[14].re);
    try std.testing.expectEqual(405, B[14].im);
    try std.testing.expectEqual(9152, B[15].re);
    try std.testing.expectEqual(432, B[15].im);
    try std.testing.expectEqual(11735, B[16].re);
    try std.testing.expectEqual(459, B[16].im);
    try std.testing.expectEqual(14318, B[17].re);
    try std.testing.expectEqual(486, B[17].im);
    try std.testing.expectEqual(513899, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(22364, B[19].re);
    try std.testing.expectEqual(540, B[19].im);
    try std.testing.expectEqual(10223, B[20].re);
    try std.testing.expectEqual(567, B[20].im);
    try std.testing.expectEqual(13346, B[21].re);
    try std.testing.expectEqual(594, B[21].im);
    try std.testing.expectEqual(16469, B[22].re);
    try std.testing.expectEqual(621, B[22].im);
    try std.testing.expectEqual(19592, B[23].re);
    try std.testing.expectEqual(648, B[23].im);
    try std.testing.expectEqual(772625, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .RowMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(252323, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(4022, B[1].re);
    try std.testing.expectEqual(54, B[1].im);
    try std.testing.expectEqual(4985, B[2].re);
    try std.testing.expectEqual(81, B[2].im);
    try std.testing.expectEqual(5948, B[3].re);
    try std.testing.expectEqual(108, B[3].im);
    try std.testing.expectEqual(6911, B[4].re);
    try std.testing.expectEqual(135, B[4].im);
    try std.testing.expectEqual(21734, B[5].re);
    try std.testing.expectEqual(486, B[5].im);
    try std.testing.expectEqual(509393, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(8096, B[7].re);
    try std.testing.expectEqual(216, B[7].im);
    try std.testing.expectEqual(10079, B[8].re);
    try std.testing.expectEqual(243, B[8].im);
    try std.testing.expectEqual(12062, B[9].re);
    try std.testing.expectEqual(270, B[9].im);
    try std.testing.expectEqual(25019, B[10].re);
    try std.testing.expectEqual(891, B[10].im);
    try std.testing.expectEqual(31244, B[11].re);
    try std.testing.expectEqual(972, B[11].im);
    try std.testing.expectEqual(939575, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(14210, B[13].re);
    try std.testing.expectEqual(378, B[13].im);
    try std.testing.expectEqual(17213, B[14].re);
    try std.testing.expectEqual(405, B[14].im);
    try std.testing.expectEqual(28304, B[15].re);
    try std.testing.expectEqual(1296, B[15].im);
    try std.testing.expectEqual(36161, B[16].re);
    try std.testing.expectEqual(1377, B[16].im);
    try std.testing.expectEqual(44018, B[17].re);
    try std.testing.expectEqual(1458, B[17].im);
    try std.testing.expectEqual(1542869, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(22364, B[19].re);
    try std.testing.expectEqual(540, B[19].im);
    try std.testing.expectEqual(31589, B[20].re);
    try std.testing.expectEqual(1701, B[20].im);
    try std.testing.expectEqual(41078, B[21].re);
    try std.testing.expectEqual(1782, B[21].im);
    try std.testing.expectEqual(50567, B[22].re);
    try std.testing.expectEqual(1863, B[22].im);
    try std.testing.expectEqual(60056, B[23].re);
    try std.testing.expectEqual(1944, B[23].im);
    try std.testing.expectEqual(2319275, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    blas.herk(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(757025, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(12194, B[1].re);
    try std.testing.expectEqual(162, B[1].im);
    try std.testing.expectEqual(15155, B[2].re);
    try std.testing.expectEqual(243, B[2].im);
    try std.testing.expectEqual(18116, B[3].re);
    try std.testing.expectEqual(324, B[3].im);
    try std.testing.expectEqual(21077, B[4].re);
    try std.testing.expectEqual(405, B[4].im);
    try std.testing.expectEqual(21734, B[5].re);
    try std.testing.expectEqual(486, B[5].im);
    try std.testing.expectEqual(1528487, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(24776, B[7].re);
    try std.testing.expectEqual(648, B[7].im);
    try std.testing.expectEqual(30905, B[8].re);
    try std.testing.expectEqual(729, B[8].im);
    try std.testing.expectEqual(37034, B[9].re);
    try std.testing.expectEqual(810, B[9].im);
    try std.testing.expectEqual(25019, B[10].re);
    try std.testing.expectEqual(891, B[10].im);
    try std.testing.expectEqual(31244, B[11].re);
    try std.testing.expectEqual(972, B[11].im);
    try std.testing.expectEqual(2819501, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(43694, B[13].re);
    try std.testing.expectEqual(1134, B[13].im);
    try std.testing.expectEqual(52991, B[14].re);
    try std.testing.expectEqual(1215, B[14].im);
    try std.testing.expectEqual(28304, B[15].re);
    try std.testing.expectEqual(1296, B[15].im);
    try std.testing.expectEqual(36161, B[16].re);
    try std.testing.expectEqual(1377, B[16].im);
    try std.testing.expectEqual(44018, B[17].re);
    try std.testing.expectEqual(1458, B[17].im);
    try std.testing.expectEqual(4630067, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(68948, B[19].re);
    try std.testing.expectEqual(1620, B[19].im);
    try std.testing.expectEqual(31589, B[20].re);
    try std.testing.expectEqual(1701, B[20].im);
    try std.testing.expectEqual(41078, B[21].re);
    try std.testing.expectEqual(1782, B[21].im);
    try std.testing.expectEqual(50567, B[22].re);
    try std.testing.expectEqual(1863, B[22].im);
    try std.testing.expectEqual(60056, B[23].re);
    try std.testing.expectEqual(1944, B[23].im);
    try std.testing.expectEqual(6960185, B[24].re);
    try std.testing.expectEqual(0, B[24].im);
}
