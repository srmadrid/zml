const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

pub inline fn syr2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or trans == .ConjTrans or trans == .ConjNoTrans) return;

    var UPLO = uplo;
    var TRANS = trans;
    var NROWAB = if (trans == .NoTrans) n else k;
    const ldcp1 = ldc + 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANS = if (trans == .NoTrans) .Trans else .NoTrans;
        NROWAB = if (trans == .NoTrans) k else n;
    }

    if (lda < @max(1, NROWAB)) return;
    if (ldb < @max(1, NROWAB)) return;
    if (ldc < @max(1, n)) return;

    switch (numericType) {
        .bool => @compileError("blas.syr2k does not support bool."),
        .int, .float => {
            if ((alpha == 0 or k <= 0) and beta == 1) return;

            if (alpha == 0) {
                if (UPLO == .Upper) {
                    if (beta == 0) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = 0;
                            var icij: isize = jcj;
                            while (i <= j) {
                                C[@intCast(icij)] = 0;

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
                            while (i <= j) {
                                C[@intCast(icij)] *= beta;

                                i += 1;
                                icij += 1;
                            }

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
                            var icij: isize = j + jcj;
                            while (i < n) {
                                C[@intCast(icij)] = 0;

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
                            var i: isize = j;
                            var icij: isize = j + jcj;
                            while (i < n) {
                                C[@intCast(icij)] *= beta;

                                i += 1;
                                icij += 1;
                            }

                            j += 1;
                            jcj += ldc;
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
                            while (icij < j + 1) {
                                Cpjcj[@intCast(icij)] = 0;

                                icij += 1;
                            }
                        } else if (beta != 1) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < j + 1) {
                                Cpjcj[@intCast(icij)] *= beta;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var ibjl: isize = ibj;
                        var jal: isize = 0;
                        var jbl: isize = 0;
                        while (l < k) {
                            const t0 = alpha * A[@intCast(iajl)];
                            const t1 = alpha * B[@intCast(ibjl)];

                            var i: isize = 0;
                            var iail: isize = jal;
                            var ibil: isize = jbl;
                            var icij: isize = jcj;
                            while (i <= j) {
                                C[@intCast(icij)] += t1 * A[@intCast(iail)] + t0 * B[@intCast(ibil)];

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
                        var i: isize = 0;
                        var jai: isize = 0;
                        var jbi: isize = 0;
                        var icij: isize = jcj;
                        while (i <= j) {
                            var t0: T = 0;
                            var t1: T = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0 += A[@intCast(iali)] * B[@intCast(iblj)];
                                t1 += B[@intCast(ibli)] * A[@intCast(ialj)];

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0 + alpha * t1;

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
                        if (beta == 0) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(j + jcj)]);
                            while (icj < n - j) {
                                Cpicij[@intCast(icj)] = 0;

                                icj += 1;
                            }
                        } else if (beta != 1) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(j + jcj)]);
                            while (icj < n - j) {
                                Cpicij[@intCast(icj)] *= beta;

                                icj += 1;
                            }
                        }

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var ibjl: isize = ibj;
                        var jal: isize = 0;
                        var jbl: isize = 0;
                        while (l < k) {
                            const t0 = alpha * A[@intCast(iajl)];
                            const t1 = alpha * B[@intCast(ibjl)];

                            var i: isize = j;
                            var iail: isize = j + jal;
                            var ibil: isize = j + jbl;
                            var icij: isize = j + jcj;
                            while (i < n) {
                                C[@intCast(icij)] += t1 * A[@intCast(iail)] + t0 * B[@intCast(ibil)];

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
                            var t0: T = 0;
                            var t1: T = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0 += A[@intCast(iali)] * B[@intCast(iblj)];
                                t1 += B[@intCast(ibli)] * A[@intCast(ialj)];

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0 + alpha * t1;

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
        .cfloat => {
            if (((alpha.re == 0 and alpha.im == 0) or k <= 0) and beta.re == 1 and beta.im == 0) return;

            if (alpha.re == 0 and alpha.im == 0) {
                if (UPLO == .Upper) {
                    if (beta.re == 0 and beta.im == 0) {
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
                    } else if (beta.re != 1 or beta.im == 0) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = 0;
                            var icij: isize = jcj;
                            while (i <= j) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;

                                i += 1;
                                icij += 1;
                            }

                            j += 1;
                            jcj += ldc;
                        }
                    }
                } else {
                    if (beta.re == 0 and beta.im == 0) {
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
                    } else if (beta.re != 1 or beta.im == 0) {
                        var j: isize = 0;
                        var jcj: isize = 0;
                        while (j < n) {
                            var i: isize = j;
                            var icij: isize = jcj;
                            while (i < n) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;

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
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < j + 1) {
                                Cpjcj[@intCast(icij)].re = 0;
                                Cpjcj[@intCast(icij)].im = 0;

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < j + 1) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var ibjl: isize = ibj;
                        var jal: isize = 0;
                        var jbl: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha.re * B[@intCast(ibjl)].re - alpha.im * B[@intCast(ibjl)].im;
                            t0.im = alpha.im * B[@intCast(ibjl)].re + alpha.re * B[@intCast(ibjl)].im;
                            var t1: T = undefined;
                            t1.re = alpha.re * A[@intCast(iajl)].re - alpha.im * A[@intCast(iajl)].im;
                            t1.im = alpha.re * A[@intCast(iajl)].im + alpha.im * A[@intCast(iajl)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var ibil: isize = jbl;
                            var icij: isize = jcj;
                            while (i <= j) {
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
                        var i: isize = 0;
                        var jai: isize = 0;
                        var jbi: isize = 0;
                        var icij: isize = jcj;
                        while (i <= j) {
                            var t0: T = undefined;
                            t0.re = 0;
                            t0.im = 0;
                            var t1: T = undefined;
                            t1.re = 0;
                            t1.im = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * B[@intCast(iblj)].re - A[@intCast(iali)].im * B[@intCast(iblj)].im;
                                t0.im += A[@intCast(iali)].re * B[@intCast(iblj)].im + A[@intCast(iali)].im * B[@intCast(iblj)].re;
                                t1.re += B[@intCast(ibli)].re * A[@intCast(ialj)].re - B[@intCast(ibli)].im * A[@intCast(ialj)].im;
                                t1.im += B[@intCast(ibli)].re * A[@intCast(ialj)].im + B[@intCast(ibli)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im == 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im + alpha.re * t1.re - alpha.im * t1.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re + alpha.re * t1.im + alpha.im * t1.re;

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
                        if (beta.re == 0 and beta.im == 0) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(j + jcj)]);
                            while (icj < n - j) {
                                Cpicij[@intCast(icj)].re = 0;
                                Cpicij[@intCast(icj)].im = 0;

                                icj += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(j + jcj)]);
                            while (icj < n - j) {
                                const tmp = Cpicij[@intCast(icj)].re * beta.re - Cpicij[@intCast(icj)].im * beta.im;
                                Cpicij[@intCast(icj)].im = Cpicij[@intCast(icj)].re * beta.im + Cpicij[@intCast(icj)].im * beta.re;
                                Cpicij[@intCast(icj)].re = tmp;

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
                            t0.re = alpha.re * B[@intCast(ibjl)].re - alpha.im * B[@intCast(ibjl)].im;
                            t0.im = alpha.im * B[@intCast(ibjl)].re + alpha.re * B[@intCast(ibjl)].im;
                            var t1: T = undefined;
                            t1.re = alpha.re * A[@intCast(iajl)].re - alpha.im * A[@intCast(iajl)].im;
                            t1.im = alpha.re * A[@intCast(iajl)].im + alpha.im * A[@intCast(iajl)].re;

                            var i: isize = j;
                            var iail: isize = j + jal;
                            var ibil: isize = j + jbl;
                            var icij: isize = j + jcj;
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
                            var t0: T = undefined;
                            t0.re = 0;
                            t0.im = 0;
                            var t1: T = undefined;
                            t1.re = 0;
                            t1.im = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            var ibli: isize = jbi;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * B[@intCast(iblj)].re - A[@intCast(iali)].im * B[@intCast(iblj)].im;
                                t0.im += A[@intCast(iali)].re * B[@intCast(iblj)].im + A[@intCast(iali)].im * B[@intCast(iblj)].re;
                                t1.re += A[@intCast(ialj)].re * B[@intCast(ibli)].re - A[@intCast(ialj)].im * B[@intCast(ibli)].im;
                                t1.im += A[@intCast(ialj)].re * B[@intCast(ibli)].im + A[@intCast(ialj)].im * B[@intCast(ibli)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                                ibli += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im + alpha.re * t1.re - alpha.im * t1.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re + alpha.re * t1.im + alpha.im * t1.re;

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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.syr2k only supports simple types."),
        .unsupported => unreachable,
    }
}

test syr2k {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 3;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(f64, n * k);
    defer a.free(A);
    const B = try a.alloc(f64, n * k);
    defer a.free(B);
    const C = try a.alloc(f64, n * n);
    defer a.free(C);

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
    });
    @memcpy(B.ptr, &[_]f64{
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
    });
    @memcpy(C.ptr, &[_]f64{
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
        21,
        22,
        23,
        24,
        25,
    });

    blas.syr2k(f64, .RowMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(59, C[0]);
    try std.testing.expectEqual(134, C[1]);
    try std.testing.expectEqual(209, C[2]);
    try std.testing.expectEqual(284, C[3]);
    try std.testing.expectEqual(359, C[4]);
    try std.testing.expectEqual(6, C[5]);
    try std.testing.expectEqual(329, C[6]);
    try std.testing.expectEqual(512, C[7]);
    try std.testing.expectEqual(695, C[8]);
    try std.testing.expectEqual(878, C[9]);
    try std.testing.expectEqual(11, C[10]);
    try std.testing.expectEqual(12, C[11]);
    try std.testing.expectEqual(815, C[12]);
    try std.testing.expectEqual(1106, C[13]);
    try std.testing.expectEqual(1397, C[14]);
    try std.testing.expectEqual(16, C[15]);
    try std.testing.expectEqual(17, C[16]);
    try std.testing.expectEqual(18, C[17]);
    try std.testing.expectEqual(1517, C[18]);
    try std.testing.expectEqual(1916, C[19]);
    try std.testing.expectEqual(21, C[20]);
    try std.testing.expectEqual(22, C[21]);
    try std.testing.expectEqual(23, C[22]);
    try std.testing.expectEqual(24, C[23]);
    try std.testing.expectEqual(2435, C[24]);

    blas.syr2k(f64, .ColumnMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(809, C[0]);
    try std.testing.expectEqual(134, C[1]);
    try std.testing.expectEqual(209, C[2]);
    try std.testing.expectEqual(284, C[3]);
    try std.testing.expectEqual(359, C[4]);
    try std.testing.expectEqual(722, C[5]);
    try std.testing.expectEqual(1775, C[6]);
    try std.testing.expectEqual(512, C[7]);
    try std.testing.expectEqual(695, C[8]);
    try std.testing.expectEqual(878, C[9]);
    try std.testing.expectEqual(809, C[10]);
    try std.testing.expectEqual(908, C[11]);
    try std.testing.expectEqual(3413, C[12]);
    try std.testing.expectEqual(1106, C[13]);
    try std.testing.expectEqual(1397, C[14]);
    try std.testing.expectEqual(896, C[15]);
    try std.testing.expectEqual(1007, C[16]);
    try std.testing.expectEqual(1118, C[17]);
    try std.testing.expectEqual(5723, C[18]);
    try std.testing.expectEqual(1916, C[19]);
    try std.testing.expectEqual(983, C[20]);
    try std.testing.expectEqual(1106, C[21]);
    try std.testing.expectEqual(1229, C[22]);
    try std.testing.expectEqual(1352, C[23]);
    try std.testing.expectEqual(8705, C[24]);

    blas.syr2k(f64, .RowMajor, .Upper, .Trans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(3059, C[0]);
    try std.testing.expectEqual(1106, C[1]);
    try std.testing.expectEqual(1403, C[2]);
    try std.testing.expectEqual(1700, C[3]);
    try std.testing.expectEqual(1997, C[4]);
    try std.testing.expectEqual(722, C[5]);
    try std.testing.expectEqual(6113, C[6]);
    try std.testing.expectEqual(2408, C[7]);
    try std.testing.expectEqual(3041, C[8]);
    try std.testing.expectEqual(3674, C[9]);
    try std.testing.expectEqual(809, C[10]);
    try std.testing.expectEqual(908, C[11]);
    try std.testing.expectEqual(11207, C[12]);
    try std.testing.expectEqual(4382, C[13]);
    try std.testing.expectEqual(5351, C[14]);
    try std.testing.expectEqual(896, C[15]);
    try std.testing.expectEqual(1007, C[16]);
    try std.testing.expectEqual(1118, C[17]);
    try std.testing.expectEqual(18341, C[18]);
    try std.testing.expectEqual(7028, C[19]);
    try std.testing.expectEqual(983, C[20]);
    try std.testing.expectEqual(1106, C[21]);
    try std.testing.expectEqual(1229, C[22]);
    try std.testing.expectEqual(1352, C[23]);
    try std.testing.expectEqual(27515, C[24]);

    blas.syr2k(f64, .ColumnMajor, .Upper, .Trans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(9233, C[0]);
    try std.testing.expectEqual(1106, C[1]);
    try std.testing.expectEqual(1403, C[2]);
    try std.testing.expectEqual(1700, C[3]);
    try std.testing.expectEqual(1997, C[4]);
    try std.testing.expectEqual(2294, C[5]);
    try std.testing.expectEqual(18647, C[6]);
    try std.testing.expectEqual(2408, C[7]);
    try std.testing.expectEqual(3041, C[8]);
    try std.testing.expectEqual(3674, C[9]);
    try std.testing.expectEqual(2627, C[10]);
    try std.testing.expectEqual(3212, C[11]);
    try std.testing.expectEqual(34397, C[12]);
    try std.testing.expectEqual(4382, C[13]);
    try std.testing.expectEqual(5351, C[14]);
    try std.testing.expectEqual(2960, C[15]);
    try std.testing.expectEqual(3689, C[16]);
    try std.testing.expectEqual(4418, C[17]);
    try std.testing.expectEqual(56483, C[18]);
    try std.testing.expectEqual(7028, C[19]);
    try std.testing.expectEqual(3293, C[20]);
    try std.testing.expectEqual(4166, C[21]);
    try std.testing.expectEqual(5039, C[22]);
    try std.testing.expectEqual(5912, C[23]);
    try std.testing.expectEqual(84905, C[24]);

    blas.syr2k(f64, .RowMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(27755, C[0]);
    try std.testing.expectEqual(1106, C[1]);
    try std.testing.expectEqual(1403, C[2]);
    try std.testing.expectEqual(1700, C[3]);
    try std.testing.expectEqual(1997, C[4]);
    try std.testing.expectEqual(7010, C[5]);
    try std.testing.expectEqual(56249, C[6]);
    try std.testing.expectEqual(2408, C[7]);
    try std.testing.expectEqual(3041, C[8]);
    try std.testing.expectEqual(3674, C[9]);
    try std.testing.expectEqual(8081, C[10]);
    try std.testing.expectEqual(10124, C[11]);
    try std.testing.expectEqual(103967, C[12]);
    try std.testing.expectEqual(4382, C[13]);
    try std.testing.expectEqual(5351, C[14]);
    try std.testing.expectEqual(9152, C[15]);
    try std.testing.expectEqual(11735, C[16]);
    try std.testing.expectEqual(14318, C[17]);
    try std.testing.expectEqual(170909, C[18]);
    try std.testing.expectEqual(7028, C[19]);
    try std.testing.expectEqual(10223, C[20]);
    try std.testing.expectEqual(13346, C[21]);
    try std.testing.expectEqual(16469, C[22]);
    try std.testing.expectEqual(19592, C[23]);
    try std.testing.expectEqual(257075, C[24]);

    blas.syr2k(f64, .ColumnMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(83897, C[0]);
    try std.testing.expectEqual(4022, C[1]);
    try std.testing.expectEqual(4985, C[2]);
    try std.testing.expectEqual(5948, C[3]);
    try std.testing.expectEqual(6911, C[4]);
    try std.testing.expectEqual(7010, C[5]);
    try std.testing.expectEqual(169535, C[6]);
    try std.testing.expectEqual(8096, C[7]);
    try std.testing.expectEqual(10079, C[8]);
    try std.testing.expectEqual(12062, C[9]);
    try std.testing.expectEqual(8081, C[10]);
    try std.testing.expectEqual(10124, C[11]);
    try std.testing.expectEqual(312869, C[12]);
    try std.testing.expectEqual(14210, C[13]);
    try std.testing.expectEqual(17213, C[14]);
    try std.testing.expectEqual(9152, C[15]);
    try std.testing.expectEqual(11735, C[16]);
    try std.testing.expectEqual(14318, C[17]);
    try std.testing.expectEqual(513899, C[18]);
    try std.testing.expectEqual(22364, C[19]);
    try std.testing.expectEqual(10223, C[20]);
    try std.testing.expectEqual(13346, C[21]);
    try std.testing.expectEqual(16469, C[22]);
    try std.testing.expectEqual(19592, C[23]);
    try std.testing.expectEqual(772625, C[24]);

    blas.syr2k(f64, .RowMajor, .Lower, .Trans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(252323, C[0]);
    try std.testing.expectEqual(4022, C[1]);
    try std.testing.expectEqual(4985, C[2]);
    try std.testing.expectEqual(5948, C[3]);
    try std.testing.expectEqual(6911, C[4]);
    try std.testing.expectEqual(21734, C[5]);
    try std.testing.expectEqual(509393, C[6]);
    try std.testing.expectEqual(8096, C[7]);
    try std.testing.expectEqual(10079, C[8]);
    try std.testing.expectEqual(12062, C[9]);
    try std.testing.expectEqual(25019, C[10]);
    try std.testing.expectEqual(31244, C[11]);
    try std.testing.expectEqual(939575, C[12]);
    try std.testing.expectEqual(14210, C[13]);
    try std.testing.expectEqual(17213, C[14]);
    try std.testing.expectEqual(28304, C[15]);
    try std.testing.expectEqual(36161, C[16]);
    try std.testing.expectEqual(44018, C[17]);
    try std.testing.expectEqual(1542869, C[18]);
    try std.testing.expectEqual(22364, C[19]);
    try std.testing.expectEqual(31589, C[20]);
    try std.testing.expectEqual(41078, C[21]);
    try std.testing.expectEqual(50567, C[22]);
    try std.testing.expectEqual(60056, C[23]);
    try std.testing.expectEqual(2319275, C[24]);

    blas.syr2k(f64, .ColumnMajor, .Lower, .Trans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(757025, C[0]);
    try std.testing.expectEqual(12194, C[1]);
    try std.testing.expectEqual(15155, C[2]);
    try std.testing.expectEqual(18116, C[3]);
    try std.testing.expectEqual(21077, C[4]);
    try std.testing.expectEqual(21734, C[5]);
    try std.testing.expectEqual(1528487, C[6]);
    try std.testing.expectEqual(24776, C[7]);
    try std.testing.expectEqual(30905, C[8]);
    try std.testing.expectEqual(37034, C[9]);
    try std.testing.expectEqual(25019, C[10]);
    try std.testing.expectEqual(31244, C[11]);
    try std.testing.expectEqual(2819501, C[12]);
    try std.testing.expectEqual(43694, C[13]);
    try std.testing.expectEqual(52991, C[14]);
    try std.testing.expectEqual(28304, C[15]);
    try std.testing.expectEqual(36161, C[16]);
    try std.testing.expectEqual(44018, C[17]);
    try std.testing.expectEqual(4630067, C[18]);
    try std.testing.expectEqual(68948, C[19]);
    try std.testing.expectEqual(31589, C[20]);
    try std.testing.expectEqual(41078, C[21]);
    try std.testing.expectEqual(50567, C[22]);
    try std.testing.expectEqual(60056, C[23]);
    try std.testing.expectEqual(6960185, C[24]);

    const gamma = Complex(f64).init(2, 2);
    const delta = Complex(f64).init(3, 3);

    const D = try a.alloc(Complex(f64), n * k);
    defer a.free(D);
    const E = try a.alloc(Complex(f64), n * k);
    defer a.free(E);
    const F = try a.alloc(Complex(f64), n * n);
    defer a.free(F);

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
    });
    @memcpy(E.ptr, &[_]Complex(f64){
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
    @memcpy(F.ptr, &[_]Complex(f64){
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

    blas.syr2k(Complex(f64), .RowMajor, .Upper, .NoTrans, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n);

    try std.testing.expectEqual(-112, F[0].re);
    try std.testing.expectEqual(118, F[0].im);
    try std.testing.expectEqual(-256, F[1].re);
    try std.testing.expectEqual(268, F[1].im);
    try std.testing.expectEqual(-400, F[2].re);
    try std.testing.expectEqual(418, F[2].im);
    try std.testing.expectEqual(-544, F[3].re);
    try std.testing.expectEqual(568, F[3].im);
    try std.testing.expectEqual(-688, F[4].re);
    try std.testing.expectEqual(718, F[4].im);
    try std.testing.expectEqual(6, F[5].re);
    try std.testing.expectEqual(6, F[5].im);
    try std.testing.expectEqual(-616, F[6].re);
    try std.testing.expectEqual(658, F[6].im);
    try std.testing.expectEqual(-976, F[7].re);
    try std.testing.expectEqual(1024, F[7].im);
    try std.testing.expectEqual(-1336, F[8].re);
    try std.testing.expectEqual(1390, F[8].im);
    try std.testing.expectEqual(-1696, F[9].re);
    try std.testing.expectEqual(1756, F[9].im);
    try std.testing.expectEqual(11, F[10].re);
    try std.testing.expectEqual(11, F[10].im);
    try std.testing.expectEqual(12, F[11].re);
    try std.testing.expectEqual(12, F[11].im);
    try std.testing.expectEqual(-1552, F[12].re);
    try std.testing.expectEqual(1630, F[12].im);
    try std.testing.expectEqual(-2128, F[13].re);
    try std.testing.expectEqual(2212, F[13].im);
    try std.testing.expectEqual(-2704, F[14].re);
    try std.testing.expectEqual(2794, F[14].im);
    try std.testing.expectEqual(16, F[15].re);
    try std.testing.expectEqual(16, F[15].im);
    try std.testing.expectEqual(17, F[16].re);
    try std.testing.expectEqual(17, F[16].im);
    try std.testing.expectEqual(18, F[17].re);
    try std.testing.expectEqual(18, F[17].im);
    try std.testing.expectEqual(-2920, F[18].re);
    try std.testing.expectEqual(3034, F[18].im);
    try std.testing.expectEqual(-3712, F[19].re);
    try std.testing.expectEqual(3832, F[19].im);
    try std.testing.expectEqual(21, F[20].re);
    try std.testing.expectEqual(21, F[20].im);
    try std.testing.expectEqual(22, F[21].re);
    try std.testing.expectEqual(22, F[21].im);
    try std.testing.expectEqual(23, F[22].re);
    try std.testing.expectEqual(23, F[22].im);
    try std.testing.expectEqual(24, F[23].re);
    try std.testing.expectEqual(24, F[23].im);
    try std.testing.expectEqual(-4720, F[24].re);
    try std.testing.expectEqual(4870, F[24].im);

    blas.syr2k(Complex(f64), .ColumnMajor, .Upper, .NoTrans, n, k, gamma, D.ptr, n, E.ptr, n, delta, F.ptr, n);

    try std.testing.expectEqual(-1954, F[0].re);
    try std.testing.expectEqual(1282, F[0].im);
    try std.testing.expectEqual(-256, F[1].re);
    try std.testing.expectEqual(268, F[1].im);
    try std.testing.expectEqual(-400, F[2].re);
    try std.testing.expectEqual(418, F[2].im);
    try std.testing.expectEqual(-544, F[3].re);
    try std.testing.expectEqual(568, F[3].im);
    try std.testing.expectEqual(-688, F[4].re);
    try std.testing.expectEqual(718, F[4].im);
    try std.testing.expectEqual(-1408, F[5].re);
    try std.testing.expectEqual(1444, F[5].im);
    try std.testing.expectEqual(-5398, F[6].re);
    try std.testing.expectEqual(1702, F[6].im);
    try std.testing.expectEqual(-976, F[7].re);
    try std.testing.expectEqual(1024, F[7].im);
    try std.testing.expectEqual(-1336, F[8].re);
    try std.testing.expectEqual(1390, F[8].im);
    try std.testing.expectEqual(-1696, F[9].re);
    try std.testing.expectEqual(1756, F[9].im);
    try std.testing.expectEqual(-1552, F[10].re);
    try std.testing.expectEqual(1618, F[10].im);
    try std.testing.expectEqual(-1744, F[11].re);
    try std.testing.expectEqual(1816, F[11].im);
    try std.testing.expectEqual(-11482, F[12].re);
    try std.testing.expectEqual(2170, F[12].im);
    try std.testing.expectEqual(-2128, F[13].re);
    try std.testing.expectEqual(2212, F[13].im);
    try std.testing.expectEqual(-2704, F[14].re);
    try std.testing.expectEqual(2794, F[14].im);
    try std.testing.expectEqual(-1696, F[15].re);
    try std.testing.expectEqual(1792, F[15].im);
    try std.testing.expectEqual(-1912, F[16].re);
    try std.testing.expectEqual(2014, F[16].im);
    try std.testing.expectEqual(-2128, F[17].re);
    try std.testing.expectEqual(2236, F[17].im);
    try std.testing.expectEqual(-20206, F[18].re);
    try std.testing.expectEqual(2686, F[18].im);
    try std.testing.expectEqual(-3712, F[19].re);
    try std.testing.expectEqual(3832, F[19].im);
    try std.testing.expectEqual(-1840, F[20].re);
    try std.testing.expectEqual(1966, F[20].im);
    try std.testing.expectEqual(-2080, F[21].re);
    try std.testing.expectEqual(2212, F[21].im);
    try std.testing.expectEqual(-2320, F[22].re);
    try std.testing.expectEqual(2458, F[22].im);
    try std.testing.expectEqual(-2560, F[23].re);
    try std.testing.expectEqual(2704, F[23].im);
    try std.testing.expectEqual(-31570, F[24].re);
    try std.testing.expectEqual(3250, F[24].im);

    blas.syr2k(Complex(f64), .RowMajor, .Upper, .Trans, n, k, gamma, D.ptr, n, E.ptr, n, delta, F.ptr, n);

    try std.testing.expectEqual(-10972, F[0].re);
    try std.testing.expectEqual(-752, F[0].im);
    try std.testing.expectEqual(-2980, F[1].re);
    try std.testing.expectEqual(1444, F[1].im);
    try std.testing.expectEqual(-4006, F[2].re);
    try std.testing.expectEqual(1606, F[2].im);
    try std.testing.expectEqual(-5032, F[3].re);
    try std.testing.expectEqual(1768, F[3].im);
    try std.testing.expectEqual(-6058, F[4].re);
    try std.testing.expectEqual(1930, F[4].im);
    try std.testing.expectEqual(-1408, F[5].re);
    try std.testing.expectEqual(1444, F[5].im);
    try std.testing.expectEqual(-22876, F[6].re);
    try std.testing.expectEqual(-9512, F[6].im);
    try std.testing.expectEqual(-7744, F[7].re);
    try std.testing.expectEqual(1888, F[7].im);
    try std.testing.expectEqual(-10090, F[8].re);
    try std.testing.expectEqual(2074, F[8].im);
    try std.testing.expectEqual(-12436, F[9].re);
    try std.testing.expectEqual(2260, F[9].im);
    try std.testing.expectEqual(-1552, F[10].re);
    try std.testing.expectEqual(1618, F[10].im);
    try std.testing.expectEqual(-1744, F[11].re);
    try std.testing.expectEqual(1816, F[11].im);
    try std.testing.expectEqual(-42892, F[12].re);
    try std.testing.expectEqual(-26000, F[12].im);
    try std.testing.expectEqual(-15148, F[13].re);
    try std.testing.expectEqual(2380, F[13].im);
    try std.testing.expectEqual(-18814, F[14].re);
    try std.testing.expectEqual(2590, F[14].im);
    try std.testing.expectEqual(-1696, F[15].re);
    try std.testing.expectEqual(1792, F[15].im);
    try std.testing.expectEqual(-1912, F[16].re);
    try std.testing.expectEqual(2014, F[16].im);
    try std.testing.expectEqual(-2128, F[17].re);
    try std.testing.expectEqual(2236, F[17].im);
    try std.testing.expectEqual(-71020, F[18].re);
    try std.testing.expectEqual(-50216, F[18].im);
    try std.testing.expectEqual(-25192, F[19].re);
    try std.testing.expectEqual(2920, F[19].im);
    try std.testing.expectEqual(-1840, F[20].re);
    try std.testing.expectEqual(1966, F[20].im);
    try std.testing.expectEqual(-2080, F[21].re);
    try std.testing.expectEqual(2212, F[21].im);
    try std.testing.expectEqual(-2320, F[22].re);
    try std.testing.expectEqual(2458, F[22].im);
    try std.testing.expectEqual(-2560, F[23].re);
    try std.testing.expectEqual(2704, F[23].im);
    try std.testing.expectEqual(-107260, F[24].re);
    try std.testing.expectEqual(-82160, F[24].im);

    blas.syr2k(Complex(f64), .ColumnMajor, .Upper, .Trans, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n);

    try std.testing.expectEqual(-30772, F[0].re);
    try std.testing.expectEqual(-35060, F[0].im);
    try std.testing.expectEqual(-2980, F[1].re);
    try std.testing.expectEqual(1444, F[1].im);
    try std.testing.expectEqual(-4006, F[2].re);
    try std.testing.expectEqual(1606, F[2].im);
    try std.testing.expectEqual(-5032, F[3].re);
    try std.testing.expectEqual(1768, F[3].im);
    try std.testing.expectEqual(-6058, F[4].re);
    try std.testing.expectEqual(1930, F[4].im);
    try std.testing.expectEqual(-8812, F[5].re);
    try std.testing.expectEqual(364, F[5].im);
    try std.testing.expectEqual(-40708, F[6].re);
    try std.testing.expectEqual(-96548, F[6].im);
    try std.testing.expectEqual(-7744, F[7].re);
    try std.testing.expectEqual(1888, F[7].im);
    try std.testing.expectEqual(-10090, F[8].re);
    try std.testing.expectEqual(2074, F[8].im);
    try std.testing.expectEqual(-12436, F[9].re);
    try std.testing.expectEqual(2260, F[9].im);
    try std.testing.expectEqual(-9910, F[10].re);
    try std.testing.expectEqual(598, F[10].im);
    try std.testing.expectEqual(-11656, F[11].re);
    try std.testing.expectEqual(1192, F[11].im);
    try std.testing.expectEqual(-52228, F[12].re);
    try std.testing.expectEqual(-205124, F[12].im);
    try std.testing.expectEqual(-15148, F[13].re);
    try std.testing.expectEqual(2380, F[13].im);
    try std.testing.expectEqual(-18814, F[14].re);
    try std.testing.expectEqual(2590, F[14].im);
    try std.testing.expectEqual(-11008, F[15].re);
    try std.testing.expectEqual(832, F[15].im);
    try std.testing.expectEqual(-13114, F[16].re);
    try std.testing.expectEqual(1642, F[16].im);
    try std.testing.expectEqual(-15220, F[17].re);
    try std.testing.expectEqual(2452, F[17].im);
    try std.testing.expectEqual(-65332, F[18].re);
    try std.testing.expectEqual(-360788, F[18].im);
    try std.testing.expectEqual(-25192, F[19].re);
    try std.testing.expectEqual(2920, F[19].im);
    try std.testing.expectEqual(-12106, F[20].re);
    try std.testing.expectEqual(1066, F[20].im);
    try std.testing.expectEqual(-14572, F[21].re);
    try std.testing.expectEqual(2092, F[21].im);
    try std.testing.expectEqual(-17038, F[22].re);
    try std.testing.expectEqual(3118, F[22].im);
    try std.testing.expectEqual(-19504, F[23].re);
    try std.testing.expectEqual(4144, F[23].im);
    try std.testing.expectEqual(-80020, F[24].re);
    try std.testing.expectEqual(-563540, F[24].im);

    blas.syr2k(Complex(f64), .RowMajor, .Lower, .NoTrans, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n);

    try std.testing.expectEqual(12752, F[0].re);
    try std.testing.expectEqual(-197384, F[0].im);
    try std.testing.expectEqual(-2980, F[1].re);
    try std.testing.expectEqual(1444, F[1].im);
    try std.testing.expectEqual(-4006, F[2].re);
    try std.testing.expectEqual(1606, F[2].im);
    try std.testing.expectEqual(-5032, F[3].re);
    try std.testing.expectEqual(1768, F[3].im);
    try std.testing.expectEqual(-6058, F[4].re);
    try std.testing.expectEqual(1930, F[4].im);
    try std.testing.expectEqual(-27784, F[5].re);
    try std.testing.expectEqual(-25088, F[5].im);
    try std.testing.expectEqual(166904, F[6].re);
    try std.testing.expectEqual(-411152, F[6].im);
    try std.testing.expectEqual(-7744, F[7].re);
    try std.testing.expectEqual(1888, F[7].im);
    try std.testing.expectEqual(-10090, F[8].re);
    try std.testing.expectEqual(2074, F[8].im);
    try std.testing.expectEqual(-12436, F[9].re);
    try std.testing.expectEqual(2260, F[9].im);
    try std.testing.expectEqual(-31924, F[10].re);
    try std.testing.expectEqual(-27536, F[10].im);
    try std.testing.expectEqual(-39520, F[11].re);
    try std.testing.expectEqual(-30416, F[11].im);
    try std.testing.expectEqual(457136, F[12].re);
    try std.testing.expectEqual(-770504, F[12].im);
    try std.testing.expectEqual(-15148, F[13].re);
    try std.testing.expectEqual(2380, F[13].im);
    try std.testing.expectEqual(-18814, F[14].re);
    try std.testing.expectEqual(2590, F[14].im);
    try std.testing.expectEqual(-36064, F[15].re);
    try std.testing.expectEqual(-29984, F[15].im);
    try std.testing.expectEqual(-45604, F[16].re);
    try std.testing.expectEqual(-33080, F[16].im);
    try std.testing.expectEqual(-55144, F[17].re);
    try std.testing.expectEqual(-36176, F[17].im);
    try std.testing.expectEqual(883448, F[18].re);
    try std.testing.expectEqual(-1275440, F[18].im);
    try std.testing.expectEqual(-25192, F[19].re);
    try std.testing.expectEqual(2920, F[19].im);
    try std.testing.expectEqual(-40204, F[20].re);
    try std.testing.expectEqual(-32432, F[20].im);
    try std.testing.expectEqual(-51688, F[21].re);
    try std.testing.expectEqual(-35744, F[21].im);
    try std.testing.expectEqual(-63172, F[22].re);
    try std.testing.expectEqual(-39056, F[22].im);
    try std.testing.expectEqual(-74656, F[23].re);
    try std.testing.expectEqual(-42368, F[23].im);
    try std.testing.expectEqual(1445840, F[24].re);
    try std.testing.expectEqual(-1925960, F[24].im);

    blas.syr2k(Complex(f64), .ColumnMajor, .Lower, .NoTrans, n, k, gamma, D.ptr, n, E.ptr, n, delta, F.ptr, n);

    try std.testing.expectEqual(629144, F[0].re);
    try std.testing.expectEqual(-552632, F[0].im);
    try std.testing.expectEqual(-14680, F[1].re);
    try std.testing.expectEqual(-3200, F[1].im);
    try std.testing.expectEqual(-18388, F[2].re);
    try std.testing.expectEqual(-5648, F[2].im);
    try std.testing.expectEqual(-22096, F[3].re);
    try std.testing.expectEqual(-8096, F[3].im);
    try std.testing.expectEqual(-25804, F[4].re);
    try std.testing.expectEqual(-10544, F[4].im);
    try std.testing.expectEqual(-27784, F[5].re);
    try std.testing.expectEqual(-25088, F[5].im);
    try std.testing.expectEqual(1732592, F[6].re);
    try std.testing.expectEqual(-731168, F[6].im);
    try std.testing.expectEqual(-30640, F[7].re);
    try std.testing.expectEqual(-15824, F[7].im);
    try std.testing.expectEqual(-38404, F[8].re);
    try std.testing.expectEqual(-22136, F[8].im);
    try std.testing.expectEqual(-46168, F[9].re);
    try std.testing.expectEqual(-28448, F[9].im);
    try std.testing.expectEqual(-31924, F[10].re);
    try std.testing.expectEqual(-27536, F[10].im);
    try std.testing.expectEqual(-39520, F[11].re);
    try std.testing.expectEqual(-30416, F[11].im);
    try std.testing.expectEqual(3680984, F[12].re);
    try std.testing.expectEqual(-938168, F[12].im);
    try std.testing.expectEqual(-54712, F[13].re);
    try std.testing.expectEqual(-36176, F[13].im);
    try std.testing.expectEqual(-66532, F[14].re);
    try std.testing.expectEqual(-46352, F[14].im);
    try std.testing.expectEqual(-36064, F[15].re);
    try std.testing.expectEqual(-29984, F[15].im);
    try std.testing.expectEqual(-45604, F[16].re);
    try std.testing.expectEqual(-33080, F[16].im);
    try std.testing.expectEqual(-55144, F[17].re);
    try std.testing.expectEqual(-36176, F[17].im);
    try std.testing.expectEqual(6474320, F[18].re);
    try std.testing.expectEqual(-1173632, F[18].im);
    try std.testing.expectEqual(-86896, F[19].re);
    try std.testing.expectEqual(-64256, F[19].im);
    try std.testing.expectEqual(-40204, F[20].re);
    try std.testing.expectEqual(-32432, F[20].im);
    try std.testing.expectEqual(-51688, F[21].re);
    try std.testing.expectEqual(-35744, F[21].im);
    try std.testing.expectEqual(-63172, F[22].re);
    try std.testing.expectEqual(-39056, F[22].im);
    try std.testing.expectEqual(-74656, F[23].re);
    try std.testing.expectEqual(-42368, F[23].im);
    try std.testing.expectEqual(10112600, F[24].re);
    try std.testing.expectEqual(-1437560, F[24].im);

    blas.syr2k(Complex(f64), .RowMajor, .Lower, .Trans, n, k, gamma, D.ptr, n, E.ptr, n, delta, F.ptr, n);

    try std.testing.expectEqual(3544064, F[0].re);
    try std.testing.expectEqual(230800, F[0].im);
    try std.testing.expectEqual(-14680, F[1].re);
    try std.testing.expectEqual(-3200, F[1].im);
    try std.testing.expectEqual(-18388, F[2].re);
    try std.testing.expectEqual(-5648, F[2].im);
    try std.testing.expectEqual(-22096, F[3].re);
    try std.testing.expectEqual(-8096, F[3].im);
    try std.testing.expectEqual(-25804, F[4].re);
    try std.testing.expectEqual(-10544, F[4].im);
    try std.testing.expectEqual(-9496, F[5].re);
    try std.testing.expectEqual(-157208, F[5].im);
    try std.testing.expectEqual(7389704, F[6].re);
    try std.testing.expectEqual(3005848, F[6].im);
    try std.testing.expectEqual(-30640, F[7].re);
    try std.testing.expectEqual(-15824, F[7].im);
    try std.testing.expectEqual(-38404, F[8].re);
    try std.testing.expectEqual(-22136, F[8].im);
    try std.testing.expectEqual(-46168, F[9].re);
    try std.testing.expectEqual(-28448, F[9].im);
    try std.testing.expectEqual(-14716, F[10].re);
    try std.testing.expectEqual(-176828, F[10].im);
    try std.testing.expectEqual(-29056, F[11].re);
    try std.testing.expectEqual(-208064, F[11].im);
    try std.testing.expectEqual(13855520, F[12].re);
    try std.testing.expectEqual(8230384, F[12].im);
    try std.testing.expectEqual(-54712, F[13].re);
    try std.testing.expectEqual(-36176, F[13].im);
    try std.testing.expectEqual(-66532, F[14].re);
    try std.testing.expectEqual(-46352, F[14].im);
    try std.testing.expectEqual(-19936, F[15].re);
    try std.testing.expectEqual(-196448, F[15].im);
    try std.testing.expectEqual(-39484, F[16].re);
    try std.testing.expectEqual(-234140, F[16].im);
    try std.testing.expectEqual(-59032, F[17].re);
    try std.testing.expectEqual(-271832, F[17].im);
    try std.testing.expectEqual(22941512, F[18].re);
    try std.testing.expectEqual(15904408, F[18].im);
    try std.testing.expectEqual(-86896, F[19].re);
    try std.testing.expectEqual(-64256, F[19].im);
    try std.testing.expectEqual(-25156, F[20].re);
    try std.testing.expectEqual(-216068, F[20].im);
    try std.testing.expectEqual(-49912, F[21].re);
    try std.testing.expectEqual(-260216, F[21].im);
    try std.testing.expectEqual(-74668, F[22].re);
    try std.testing.expectEqual(-304364, F[22].im);
    try std.testing.expectEqual(-99424, F[23].re);
    try std.testing.expectEqual(-348512, F[23].im);
    try std.testing.expectEqual(34647680, F[24].re);
    try std.testing.expectEqual(26027920, F[24].im);

    blas.syr2k(Complex(f64), .ColumnMajor, .Lower, .Trans, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n);

    try std.testing.expectEqual(9939680, F[0].re);
    try std.testing.expectEqual(11324704, F[0].im);
    try std.testing.expectEqual(-34696, F[1].re);
    try std.testing.expectEqual(-53384, F[1].im);
    try std.testing.expectEqual(-38620, F[2].re);
    try std.testing.expectEqual(-71708, F[2].im);
    try std.testing.expectEqual(-42544, F[3].re);
    try std.testing.expectEqual(-90032, F[3].im);
    try std.testing.expectEqual(-46468, F[4].re);
    try std.testing.expectEqual(-108356, F[4].im);
    try std.testing.expectEqual(-9496, F[5].re);
    try std.testing.expectEqual(-157208, F[5].im);
    try std.testing.expectEqual(13150952, F[6].re);
    try std.testing.expectEqual(31187272, F[6].im);
    try std.testing.expectEqual(-45424, F[7].re);
    try std.testing.expectEqual(-138416, F[7].im);
    try std.testing.expectEqual(-50140, F[8].re);
    try std.testing.expectEqual(-180284, F[8].im);
    try std.testing.expectEqual(-54856, F[9].re);
    try std.testing.expectEqual(-222152, F[9].im);
    try std.testing.expectEqual(-14716, F[10].re);
    try std.testing.expectEqual(-176828, F[10].im);
    try std.testing.expectEqual(-29056, F[11].re);
    try std.testing.expectEqual(-208064, F[11].im);
    try std.testing.expectEqual(16873856, F[12].re);
    try std.testing.expectEqual(66259264, F[12].im);
    try std.testing.expectEqual(-57736, F[13].re);
    try std.testing.expectEqual(-270536, F[13].im);
    try std.testing.expectEqual(-63244, F[14].re);
    try std.testing.expectEqual(-335948, F[14].im);
    try std.testing.expectEqual(-19936, F[15].re);
    try std.testing.expectEqual(-196448, F[15].im);
    try std.testing.expectEqual(-39484, F[16].re);
    try std.testing.expectEqual(-234140, F[16].im);
    try std.testing.expectEqual(-59032, F[17].re);
    try std.testing.expectEqual(-271832, F[17].im);
    try std.testing.expectEqual(21108392, F[18].re);
    try std.testing.expectEqual(116540680, F[18].im);
    try std.testing.expectEqual(-71632, F[19].re);
    try std.testing.expectEqual(-449744, F[19].im);
    try std.testing.expectEqual(-25156, F[20].re);
    try std.testing.expectEqual(-216068, F[20].im);
    try std.testing.expectEqual(-49912, F[21].re);
    try std.testing.expectEqual(-260216, F[21].im);
    try std.testing.expectEqual(-74668, F[22].re);
    try std.testing.expectEqual(-304364, F[22].im);
    try std.testing.expectEqual(-99424, F[23].re);
    try std.testing.expectEqual(-348512, F[23].im);
    try std.testing.expectEqual(25854560, F[24].re);
    try std.testing.expectEqual(182031520, F[24].im);
}
