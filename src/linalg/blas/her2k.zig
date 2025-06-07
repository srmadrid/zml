const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

const Scalar = types.Scalar;

pub inline fn her2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

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
    }
}
