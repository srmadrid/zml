const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

const Numeric = core.types.Numeric;

pub inline fn her2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: Numeric(T), C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or trans == .Trans or trans == .ConjNoTrans) return;

    var UPLO = uplo;
    var TRANS = trans;
    var NROWAB = if (trans == .NoTrans) n else k;
    const ldcp1 = ldc + 1;
    var ALPHA = alpha;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANS = if (trans == .NoTrans) .ConjTrans else .NoTrans;
        NROWAB = if (trans == .NoTrans) k else n;
        ALPHA.im = -ALPHA.im;
    }

    if (lda < @max(1, NROWAB)) return;
    if (ldb < @max(1, NROWAB)) return;
    if (ldc < @max(1, n)) return;

    switch (numericType) {
        .bool => @compileError("blas.her2k does not support bool."),
        .int, .float => @compileError("blas.her2k does not support int or float."),
        .cfloat => {
            if (((ALPHA.re == 0 and ALPHA.im == 0) or k <= 0) and beta == 1) return;

            if (ALPHA.re == 0 and ALPHA.im == 0) {
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
                    var ibj: isize = 0;
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
                        var ibjl: isize = ibj;
                        var jal: isize = 0;
                        var jbl: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = ALPHA.re * B[@intCast(ibjl)].re + ALPHA.im * B[@intCast(ibjl)].im;
                            t0.im = ALPHA.im * B[@intCast(ibjl)].re - ALPHA.re * B[@intCast(ibjl)].im;
                            var t1: T = undefined;
                            t1.re = ALPHA.re * A[@intCast(iajl)].re - ALPHA.im * A[@intCast(iajl)].im;
                            t1.im = -(ALPHA.re * A[@intCast(iajl)].im + ALPHA.im * A[@intCast(iajl)].re);

                            var i: isize = 0;
                            var iail: isize = jal;
                            var ibil: isize = jbl;
                            icij = jcj;
                            while (i < j) {
                                C[@intCast(icij)].re += t0.re * A[@intCast(iail)].re - t0.im * A[@intCast(iail)].im + t1.re * B[@intCast(ibil)].re - t1.im * B[@intCast(ibil)].im;
                                C[@intCast(icij)].im += t0.re * A[@intCast(iail)].im + t0.im * A[@intCast(iail)].re + t1.re * B[@intCast(ibil)].im + t1.im * B[@intCast(ibil)].re;

                                i += 1;
                                iail += 1;
                                ibil += 1;
                                icij += 1;
                            }

                            C[@intCast(icij)].re += t0.re * A[@intCast(iail)].re - t0.im * A[@intCast(iail)].im + t1.re * B[@intCast(ibil)].re - t1.im * B[@intCast(ibil)].im;
                            C[@intCast(icij)].im = 0;

                            l += 1;
                            iajl += lda;
                            ibjl += ldb;
                            jal += lda;
                            jbl += ldb;
                        }

                        j += 1;
                        iaj += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jaj: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        var i: isize = 0;
                        var jai: isize = 0;
                        var jbi: isize = 0;
                        var icij: isize = jcj;
                        while (i <= j) {
                            var t0 = T.init(0, 0);
                            var t1 = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * B[@intCast(iblj)].re + A[@intCast(iali)].im * B[@intCast(iblj)].im;
                                t0.im += A[@intCast(iali)].re * B[@intCast(iblj)].im - A[@intCast(iali)].im * B[@intCast(iblj)].re;
                                t1.re += B[@intCast(ibli)].re * A[@intCast(ialj)].re + B[@intCast(ibli)].im * A[@intCast(ialj)].im;
                                t1.im += B[@intCast(ibli)].re * A[@intCast(ialj)].im - B[@intCast(ibli)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (i == j) {
                                if (beta == 0) {
                                    C[@intCast(icij)].re = 0;
                                } else if (beta != 1) {
                                    C[@intCast(icij)].re *= beta;
                                }

                                C[@intCast(icij)].re += ALPHA.re * t0.re - ALPHA.im * t0.im + ALPHA.re * t1.re + ALPHA.im * t1.im;
                                C[@intCast(icij)].im = 0;
                            } else {
                                if (beta == 0) {
                                    C[@intCast(icij)].re = 0;
                                    C[@intCast(icij)].im = 0;
                                } else if (beta != 1) {
                                    C[@intCast(icij)].re *= beta;
                                    C[@intCast(icij)].im *= beta;
                                }

                                C[@intCast(icij)].re += ALPHA.re * t0.re - ALPHA.im * t0.im + ALPHA.re * t1.re + ALPHA.im * t1.im;
                                C[@intCast(icij)].im += ALPHA.re * t0.im + ALPHA.im * t0.re + ALPHA.re * t1.im - ALPHA.im * t1.re;
                            }

                            i += 1;
                            jai += lda;
                            jbi += ldb;
                            icij += 1;
                        }

                        j += 1;
                        jaj += lda;
                        jbj += ldb;
                        jcj += ldc;
                    }
                }
            } else {
                if (TRANS == .NoTrans) {
                    var j: isize = 0;
                    var iaj: isize = 0;
                    var ibj: isize = 0;
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
                        var ibjl: isize = ibj;
                        var jal: isize = 0;
                        var jbl: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = ALPHA.re * B[@intCast(ibjl)].re + ALPHA.im * B[@intCast(ibjl)].im;
                            t0.im = ALPHA.im * B[@intCast(ibjl)].re - ALPHA.re * B[@intCast(ibjl)].im;
                            var t1: T = undefined;
                            t1.re = ALPHA.re * A[@intCast(iajl)].re - ALPHA.im * A[@intCast(iajl)].im;
                            t1.im = -(ALPHA.re * A[@intCast(iajl)].im + ALPHA.im * A[@intCast(iajl)].re);

                            var iail: isize = j + jal;
                            var ibil: isize = j + jbl;
                            icij = j + jcj;

                            C[@intCast(icij)].re += t0.re * A[@intCast(iajl)].re - t0.im * A[@intCast(iajl)].im + t1.re * B[@intCast(ibjl)].re - t1.im * B[@intCast(ibjl)].im;
                            C[@intCast(icij)].im = 0;

                            iail += 1;
                            ibil += 1;
                            icij += 1;

                            var i: isize = j + 1;
                            while (i < n) {
                                C[@intCast(icij)].re += t0.re * A[@intCast(iail)].re - t0.im * A[@intCast(iail)].im + t1.re * B[@intCast(ibil)].re - t1.im * B[@intCast(ibil)].im;
                                C[@intCast(icij)].im += t0.re * A[@intCast(iail)].im + t0.im * A[@intCast(iail)].re + t1.re * B[@intCast(ibil)].im + t1.im * B[@intCast(ibil)].re;

                                i += 1;
                                iail += 1;
                                ibil += 1;
                                icij += 1;
                            }

                            l += 1;
                            iajl += lda;
                            ibjl += ldb;
                            jal += lda;
                            jbl += ldb;
                        }

                        j += 1;
                        iaj += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jaj: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < n) {
                        var i: isize = j;
                        var jai: isize = j * lda;
                        var jbi: isize = j * ldb;
                        var icij: isize = j + jcj;
                        while (i < n) {
                            var t0 = T.init(0, 0);
                            var t1 = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * B[@intCast(iblj)].re + A[@intCast(iali)].im * B[@intCast(iblj)].im;
                                t0.im += A[@intCast(iali)].re * B[@intCast(iblj)].im - A[@intCast(iali)].im * B[@intCast(iblj)].re;
                                t1.re += B[@intCast(ibli)].re * A[@intCast(ialj)].re + B[@intCast(ibli)].im * A[@intCast(ialj)].im;
                                t1.im += B[@intCast(ibli)].re * A[@intCast(ialj)].im - B[@intCast(ibli)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (i == j) {
                                if (beta == 0) {
                                    C[@intCast(icij)].re = 0;
                                } else if (beta != 1) {
                                    C[@intCast(icij)].re *= beta;
                                }

                                C[@intCast(icij)].re += ALPHA.re * t0.re - ALPHA.im * t0.im + ALPHA.re * t1.re + ALPHA.im * t1.im;
                                C[@intCast(icij)].im = 0;
                            } else {
                                if (beta == 0) {
                                    C[@intCast(icij)].re = 0;
                                    C[@intCast(icij)].im = 0;
                                } else if (beta != 1) {
                                    C[@intCast(icij)].re *= beta;
                                    C[@intCast(icij)].im *= beta;
                                }

                                C[@intCast(icij)].re += ALPHA.re * t0.re - ALPHA.im * t0.im + ALPHA.re * t1.re + ALPHA.im * t1.im;
                                C[@intCast(icij)].im += ALPHA.re * t0.im + ALPHA.im * t0.re + ALPHA.re * t1.im - ALPHA.im * t1.re;
                            }

                            i += 1;
                            jai += lda;
                            jbi += ldb;
                            icij += 1;
                        }

                        j += 1;
                        jaj += lda;
                        jbj += ldb;
                        jcj += ldc;
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.her2k only supports simple types."),
        .unsupported => unreachable,
    }
}

test her2k {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 3;
    const alpha = Complex(f64).init(2, 2);
    const beta = 3;

    const A = try a.alloc(Complex(f64), n * k);
    defer a.free(A);
    const B = try a.alloc(Complex(f64), n * k);
    defer a.free(B);
    const C = try a.alloc(Complex(f64), n * n);
    defer a.free(C);

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
        Complex(f64).init(21, 21),
        Complex(f64).init(22, 22),
        Complex(f64).init(23, 23),
        Complex(f64).init(24, 24),
        Complex(f64).init(25, 25),
    });

    blas.her2k(Complex(f64), .RowMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(115, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(262, C[1].re);
    try std.testing.expectEqual(6, C[1].im);
    try std.testing.expectEqual(409, C[2].re);
    try std.testing.expectEqual(9, C[2].im);
    try std.testing.expectEqual(556, C[3].re);
    try std.testing.expectEqual(12, C[3].im);
    try std.testing.expectEqual(703, C[4].re);
    try std.testing.expectEqual(15, C[4].im);
    try std.testing.expectEqual(6, C[5].re);
    try std.testing.expectEqual(6, C[5].im);
    try std.testing.expectEqual(637, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(1000, C[7].re);
    try std.testing.expectEqual(24, C[7].im);
    try std.testing.expectEqual(1363, C[8].re);
    try std.testing.expectEqual(27, C[8].im);
    try std.testing.expectEqual(1726, C[9].re);
    try std.testing.expectEqual(30, C[9].im);
    try std.testing.expectEqual(11, C[10].re);
    try std.testing.expectEqual(11, C[10].im);
    try std.testing.expectEqual(12, C[11].re);
    try std.testing.expectEqual(12, C[11].im);
    try std.testing.expectEqual(1591, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(2170, C[13].re);
    try std.testing.expectEqual(42, C[13].im);
    try std.testing.expectEqual(2749, C[14].re);
    try std.testing.expectEqual(45, C[14].im);
    try std.testing.expectEqual(16, C[15].re);
    try std.testing.expectEqual(16, C[15].im);
    try std.testing.expectEqual(17, C[16].re);
    try std.testing.expectEqual(17, C[16].im);
    try std.testing.expectEqual(18, C[17].re);
    try std.testing.expectEqual(18, C[17].im);
    try std.testing.expectEqual(2977, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(3772, C[19].re);
    try std.testing.expectEqual(60, C[19].im);
    try std.testing.expectEqual(21, C[20].re);
    try std.testing.expectEqual(21, C[20].im);
    try std.testing.expectEqual(22, C[21].re);
    try std.testing.expectEqual(22, C[21].im);
    try std.testing.expectEqual(23, C[22].re);
    try std.testing.expectEqual(23, C[22].im);
    try std.testing.expectEqual(24, C[23].re);
    try std.testing.expectEqual(24, C[23].im);
    try std.testing.expectEqual(4795, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .ColumnMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(1609, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(262, C[1].re);
    try std.testing.expectEqual(6, C[1].im);
    try std.testing.expectEqual(409, C[2].re);
    try std.testing.expectEqual(9, C[2].im);
    try std.testing.expectEqual(556, C[3].re);
    try std.testing.expectEqual(12, C[3].im);
    try std.testing.expectEqual(703, C[4].re);
    try std.testing.expectEqual(15, C[4].im);
    try std.testing.expectEqual(1426, C[5].re);
    try std.testing.expectEqual(18, C[5].im);
    try std.testing.expectEqual(3487, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(1000, C[7].re);
    try std.testing.expectEqual(24, C[7].im);
    try std.testing.expectEqual(1363, C[8].re);
    try std.testing.expectEqual(27, C[8].im);
    try std.testing.expectEqual(1726, C[9].re);
    try std.testing.expectEqual(30, C[9].im);
    try std.testing.expectEqual(1585, C[10].re);
    try std.testing.expectEqual(33, C[10].im);
    try std.testing.expectEqual(1780, C[11].re);
    try std.testing.expectEqual(36, C[11].im);
    try std.testing.expectEqual(6709, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(2170, C[13].re);
    try std.testing.expectEqual(42, C[13].im);
    try std.testing.expectEqual(2749, C[14].re);
    try std.testing.expectEqual(45, C[14].im);
    try std.testing.expectEqual(1744, C[15].re);
    try std.testing.expectEqual(48, C[15].im);
    try std.testing.expectEqual(1963, C[16].re);
    try std.testing.expectEqual(51, C[16].im);
    try std.testing.expectEqual(2182, C[17].re);
    try std.testing.expectEqual(54, C[17].im);
    try std.testing.expectEqual(11275, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(3772, C[19].re);
    try std.testing.expectEqual(60, C[19].im);
    try std.testing.expectEqual(1903, C[20].re);
    try std.testing.expectEqual(63, C[20].im);
    try std.testing.expectEqual(2146, C[21].re);
    try std.testing.expectEqual(66, C[21].im);
    try std.testing.expectEqual(2389, C[22].re);
    try std.testing.expectEqual(69, C[22].im);
    try std.testing.expectEqual(2632, C[23].re);
    try std.testing.expectEqual(72, C[23].im);
    try std.testing.expectEqual(17185, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .RowMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(6091, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(1426, C[5].re);
    try std.testing.expectEqual(18, C[5].im);
    try std.testing.expectEqual(12037, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(1585, C[10].re);
    try std.testing.expectEqual(33, C[10].im);
    try std.testing.expectEqual(1780, C[11].re);
    try std.testing.expectEqual(36, C[11].im);
    try std.testing.expectEqual(22063, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(1744, C[15].re);
    try std.testing.expectEqual(48, C[15].im);
    try std.testing.expectEqual(1963, C[16].re);
    try std.testing.expectEqual(51, C[16].im);
    try std.testing.expectEqual(2182, C[17].re);
    try std.testing.expectEqual(54, C[17].im);
    try std.testing.expectEqual(36169, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(1903, C[20].re);
    try std.testing.expectEqual(63, C[20].im);
    try std.testing.expectEqual(2146, C[21].re);
    try std.testing.expectEqual(66, C[21].im);
    try std.testing.expectEqual(2389, C[22].re);
    try std.testing.expectEqual(69, C[22].im);
    try std.testing.expectEqual(2632, C[23].re);
    try std.testing.expectEqual(72, C[23].im);
    try std.testing.expectEqual(54355, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(18385, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(4534, C[5].re);
    try std.testing.expectEqual(54, C[5].im);
    try std.testing.expectEqual(36727, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(5155, C[10].re);
    try std.testing.expectEqual(99, C[10].im);
    try std.testing.expectEqual(6316, C[11].re);
    try std.testing.expectEqual(108, C[11].im);
    try std.testing.expectEqual(67741, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(5776, C[15].re);
    try std.testing.expectEqual(144, C[15].im);
    try std.testing.expectEqual(7225, C[16].re);
    try std.testing.expectEqual(153, C[16].im);
    try std.testing.expectEqual(8674, C[17].re);
    try std.testing.expectEqual(162, C[17].im);
    try std.testing.expectEqual(111427, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(6397, C[20].re);
    try std.testing.expectEqual(189, C[20].im);
    try std.testing.expectEqual(8134, C[21].re);
    try std.testing.expectEqual(198, C[21].im);
    try std.testing.expectEqual(9871, C[22].re);
    try std.testing.expectEqual(207, C[22].im);
    try std.testing.expectEqual(11608, C[23].re);
    try std.testing.expectEqual(216, C[23].im);
    try std.testing.expectEqual(167785, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .RowMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(55267, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(13858, C[5].re);
    try std.testing.expectEqual(162, C[5].im);
    try std.testing.expectEqual(110797, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(15865, C[10].re);
    try std.testing.expectEqual(297, C[10].im);
    try std.testing.expectEqual(19924, C[11].re);
    try std.testing.expectEqual(324, C[11].im);
    try std.testing.expectEqual(204775, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(17872, C[15].re);
    try std.testing.expectEqual(432, C[15].im);
    try std.testing.expectEqual(23011, C[16].re);
    try std.testing.expectEqual(459, C[16].im);
    try std.testing.expectEqual(28150, C[17].re);
    try std.testing.expectEqual(486, C[17].im);
    try std.testing.expectEqual(337201, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(19879, C[20].re);
    try std.testing.expectEqual(567, C[20].im);
    try std.testing.expectEqual(26098, C[21].re);
    try std.testing.expectEqual(594, C[21].im);
    try std.testing.expectEqual(32317, C[22].re);
    try std.testing.expectEqual(621, C[22].im);
    try std.testing.expectEqual(38536, C[23].re);
    try std.testing.expectEqual(648, C[23].im);
    try std.testing.expectEqual(508075, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .ColumnMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(167065, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(7990, C[1].re);
    try std.testing.expectEqual(54, C[1].im);
    try std.testing.expectEqual(9889, C[2].re);
    try std.testing.expectEqual(81, C[2].im);
    try std.testing.expectEqual(11788, C[3].re);
    try std.testing.expectEqual(108, C[3].im);
    try std.testing.expectEqual(13687, C[4].re);
    try std.testing.expectEqual(135, C[4].im);
    try std.testing.expectEqual(13858, C[5].re);
    try std.testing.expectEqual(162, C[5].im);
    try std.testing.expectEqual(333967, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(15976, C[7].re);
    try std.testing.expectEqual(216, C[7].im);
    try std.testing.expectEqual(19915, C[8].re);
    try std.testing.expectEqual(243, C[8].im);
    try std.testing.expectEqual(23854, C[9].re);
    try std.testing.expectEqual(270, C[9].im);
    try std.testing.expectEqual(15865, C[10].re);
    try std.testing.expectEqual(297, C[10].im);
    try std.testing.expectEqual(19924, C[11].re);
    try std.testing.expectEqual(324, C[11].im);
    try std.testing.expectEqual(616261, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(28042, C[13].re);
    try std.testing.expectEqual(378, C[13].im);
    try std.testing.expectEqual(34021, C[14].re);
    try std.testing.expectEqual(405, C[14].im);
    try std.testing.expectEqual(17872, C[15].re);
    try std.testing.expectEqual(432, C[15].im);
    try std.testing.expectEqual(23011, C[16].re);
    try std.testing.expectEqual(459, C[16].im);
    try std.testing.expectEqual(28150, C[17].re);
    try std.testing.expectEqual(486, C[17].im);
    try std.testing.expectEqual(1013947, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(44188, C[19].re);
    try std.testing.expectEqual(540, C[19].im);
    try std.testing.expectEqual(19879, C[20].re);
    try std.testing.expectEqual(567, C[20].im);
    try std.testing.expectEqual(26098, C[21].re);
    try std.testing.expectEqual(594, C[21].im);
    try std.testing.expectEqual(32317, C[22].re);
    try std.testing.expectEqual(621, C[22].im);
    try std.testing.expectEqual(38536, C[23].re);
    try std.testing.expectEqual(648, C[23].im);
    try std.testing.expectEqual(1527025, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .RowMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(502459, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(7990, C[1].re);
    try std.testing.expectEqual(54, C[1].im);
    try std.testing.expectEqual(9889, C[2].re);
    try std.testing.expectEqual(81, C[2].im);
    try std.testing.expectEqual(11788, C[3].re);
    try std.testing.expectEqual(108, C[3].im);
    try std.testing.expectEqual(13687, C[4].re);
    try std.testing.expectEqual(135, C[4].im);
    try std.testing.expectEqual(42982, C[5].re);
    try std.testing.expectEqual(486, C[5].im);
    try std.testing.expectEqual(1003477, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(15976, C[7].re);
    try std.testing.expectEqual(216, C[7].im);
    try std.testing.expectEqual(19915, C[8].re);
    try std.testing.expectEqual(243, C[8].im);
    try std.testing.expectEqual(23854, C[9].re);
    try std.testing.expectEqual(270, C[9].im);
    try std.testing.expectEqual(49147, C[10].re);
    try std.testing.expectEqual(891, C[10].im);
    try std.testing.expectEqual(61516, C[11].re);
    try std.testing.expectEqual(972, C[11].im);
    try std.testing.expectEqual(1850719, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(28042, C[13].re);
    try std.testing.expectEqual(378, C[13].im);
    try std.testing.expectEqual(34021, C[14].re);
    try std.testing.expectEqual(405, C[14].im);
    try std.testing.expectEqual(55312, C[15].re);
    try std.testing.expectEqual(1296, C[15].im);
    try std.testing.expectEqual(70945, C[16].re);
    try std.testing.expectEqual(1377, C[16].im);
    try std.testing.expectEqual(86578, C[17].re);
    try std.testing.expectEqual(1458, C[17].im);
    try std.testing.expectEqual(3044185, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(44188, C[19].re);
    try std.testing.expectEqual(540, C[19].im);
    try std.testing.expectEqual(61477, C[20].re);
    try std.testing.expectEqual(1701, C[20].im);
    try std.testing.expectEqual(80374, C[21].re);
    try std.testing.expectEqual(1782, C[21].im);
    try std.testing.expectEqual(99271, C[22].re);
    try std.testing.expectEqual(1863, C[22].im);
    try std.testing.expectEqual(118168, C[23].re);
    try std.testing.expectEqual(1944, C[23].im);
    try std.testing.expectEqual(4583875, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    blas.her2k(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(1507489, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(24226, C[1].re);
    try std.testing.expectEqual(162, C[1].im);
    try std.testing.expectEqual(30067, C[2].re);
    try std.testing.expectEqual(243, C[2].im);
    try std.testing.expectEqual(35908, C[3].re);
    try std.testing.expectEqual(324, C[3].im);
    try std.testing.expectEqual(41749, C[4].re);
    try std.testing.expectEqual(405, C[4].im);
    try std.testing.expectEqual(42982, C[5].re);
    try std.testing.expectEqual(486, C[5].im);
    try std.testing.expectEqual(3011047, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(48904, C[7].re);
    try std.testing.expectEqual(648, C[7].im);
    try std.testing.expectEqual(61081, C[8].re);
    try std.testing.expectEqual(729, C[8].im);
    try std.testing.expectEqual(73258, C[9].re);
    try std.testing.expectEqual(810, C[9].im);
    try std.testing.expectEqual(49147, C[10].re);
    try std.testing.expectEqual(891, C[10].im);
    try std.testing.expectEqual(61516, C[11].re);
    try std.testing.expectEqual(972, C[11].im);
    try std.testing.expectEqual(5553709, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(86254, C[13].re);
    try std.testing.expectEqual(1134, C[13].im);
    try std.testing.expectEqual(104767, C[14].re);
    try std.testing.expectEqual(1215, C[14].im);
    try std.testing.expectEqual(55312, C[15].re);
    try std.testing.expectEqual(1296, C[15].im);
    try std.testing.expectEqual(70945, C[16].re);
    try std.testing.expectEqual(1377, C[16].im);
    try std.testing.expectEqual(86578, C[17].re);
    try std.testing.expectEqual(1458, C[17].im);
    try std.testing.expectEqual(9135475, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(136276, C[19].re);
    try std.testing.expectEqual(1620, C[19].im);
    try std.testing.expectEqual(61477, C[20].re);
    try std.testing.expectEqual(1701, C[20].im);
    try std.testing.expectEqual(80374, C[21].re);
    try std.testing.expectEqual(1782, C[21].im);
    try std.testing.expectEqual(99271, C[22].re);
    try std.testing.expectEqual(1863, C[22].im);
    try std.testing.expectEqual(118168, C[23].re);
    try std.testing.expectEqual(1944, C[23].im);
    try std.testing.expectEqual(13756345, C[24].re);
    try std.testing.expectEqual(0, C[24].im);
}
