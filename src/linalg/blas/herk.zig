const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

const Scalar = types.Scalar;

pub inline fn herk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Scalar(T), A: [*]const T, lda: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

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
    }
}
