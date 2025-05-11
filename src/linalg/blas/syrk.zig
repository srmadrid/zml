const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;

pub inline fn syrk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or trans == .ConjTrans or trans == .ConjNoTrans) return;

    var UPLO = uplo;
    var TRANS = trans;
    var NROWA = if (trans == .NoTrans) n else k;
    const ldcp1 = ldc + 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANS = if (trans == .NoTrans) .Trans else .NoTrans;
        NROWA = if (trans == .NoTrans) k else n;
    }

    if (lda < @max(1, NROWA)) return;
    if (ldc < @max(1, n)) return;

    switch (numericType) {
        .bool => @compileError("blas.syrk does not support bool."),
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
                        var jal: isize = 0;
                        while (l < k) {
                            const t0 = alpha * A[@intCast(iajl)];

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i <= j) {
                                C[@intCast(icij)] += t0 * A[@intCast(iail)];

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
                        var i: isize = 0;
                        var jai: isize = 0;
                        var icij: isize = jcj;
                        while (i <= j) {
                            var t0: T = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            while (l < k) {
                                t0 += A[@intCast(iali)] * A[@intCast(ialj)];

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0;

                            i += 1;
                            jai += lda;
                            icij += 1;
                        }

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
                        var jal: isize = 0;
                        while (l < k) {
                            const t0 = alpha * A[@intCast(iajl)];

                            var i: isize = j;
                            var iail: isize = j + jal;
                            var icij: isize = j + jcj;
                            while (i < n) {
                                C[@intCast(icij)] += t0 * A[@intCast(iail)];

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
                        var i: isize = j;
                        var jai: isize = j * lda;
                        var icij: isize = j + jcj;
                        while (i < n) {
                            var t0: T = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            while (l < k) {
                                t0 += A[@intCast(iali)] * A[@intCast(ialj)];

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0;

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
                    } else if (beta.re != 1 or beta.im != 0) {
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
                    } else if (beta.re != 1 or beta.im != 0) {
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
                        var jal: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha.re * A[@intCast(iajl)].re - alpha.im * A[@intCast(iajl)].im;
                            t0.im = alpha.re * A[@intCast(iajl)].im + alpha.im * A[@intCast(iajl)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i <= j) {
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
                        var i: isize = 0;
                        var jai: isize = 0;
                        var icij: isize = jcj;
                        while (i <= j) {
                            var t0: T = undefined;
                            t0.re = 0;
                            t0.im = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re - A[@intCast(iali)].im * A[@intCast(ialj)].im;
                                t0.im += A[@intCast(iali)].re * A[@intCast(ialj)].im + A[@intCast(iali)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            jai += lda;
                            icij += 1;
                        }

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

                        if (beta.re == 0 and beta.im == 0) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(icij)]);
                            while (icj < n - j) {
                                Cpicij[@intCast(icj)].re = 0;
                                Cpicij[@intCast(icj)].im = 0;

                                icj += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icj: isize = 0;
                            const Cpicij: [*]T = @ptrCast(&C[@intCast(icij)]);
                            while (icj < n - j) {
                                const tmp = Cpicij[@intCast(icj)].re * beta.re - Cpicij[@intCast(icj)].im * beta.im;
                                Cpicij[@intCast(icj)].im = Cpicij[@intCast(icj)].re * beta.im + Cpicij[@intCast(icj)].im * beta.re;
                                Cpicij[@intCast(icj)].re = tmp;

                                icj += 1;
                            }
                        }

                        var l: isize = 0;
                        var iajl: isize = iaj;
                        var jal: isize = 0;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha.re * A[@intCast(iajl)].re - alpha.im * A[@intCast(iajl)].im;
                            t0.im = alpha.re * A[@intCast(iajl)].im + alpha.im * A[@intCast(iajl)].re;

                            var i: isize = j;
                            var iail: isize = j + jal;
                            icij = j + jcj;
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
                        var i: isize = j;
                        var jai: isize = j * lda;
                        var icij: isize = j + jcj;
                        while (i < n) {
                            var t0: T = undefined;
                            t0.re = 0;
                            t0.im = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ialj: isize = jaj;
                            while (l < k) {
                                t0.re += A[@intCast(iali)].re * A[@intCast(ialj)].re - A[@intCast(iali)].im * A[@intCast(ialj)].im;
                                t0.im += A[@intCast(iali)].re * A[@intCast(ialj)].im + A[@intCast(iali)].im * A[@intCast(ialj)].re;

                                l += 1;
                                iali += 1;
                                ialj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)].re = 0;
                                C[@intCast(icij)].im = 0;
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.syrk only supports simple types."),
        .unsupported => unreachable,
    }
}
