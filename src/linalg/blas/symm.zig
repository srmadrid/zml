const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Side = blas.Side;
const Uplo = blas.Uplo;

pub inline fn symm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

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

    switch (numericType) {
        .bool => @compileError("blas.symm does not support bool."),
        .int, .float => {
            if (alpha == 0 and beta == 1) return;

            if (alpha == 0) {
                if (beta == 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)] = 0;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                } else if (beta != 1) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)] *= beta;

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
                            const t0 = alpha * B[@intCast(ibij)];
                            var t1: T = 0;

                            var k: isize = 0;
                            var iaki: isize = jai;
                            var ibkj: isize = jbj;
                            var ickj: isize = jcj;
                            while (k < i) {
                                C[@intCast(ickj)] += t0 * A[@intCast(iaki)];
                                t1 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                k += 1;
                                iaki += 1;
                                ibkj += 1;
                                ickj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += t0 * A[@intCast(i + jai)] + t1 * alpha;

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
                            const t0 = alpha * B[@intCast(ibij)];
                            var t1: T = 0;

                            var k: isize = i + 1;
                            var iaki: isize = jai + i + 1;
                            var ibkj: isize = jbj + i + 1;
                            var ickj: isize = jcj + i + 1;
                            while (k < M) {
                                C[@intCast(ickj)] += t0 * A[@intCast(iaki)];
                                t1 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                k += 1;
                                iaki += 1;
                                ibkj += 1;
                                ickj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += t0 * A[@intCast(i + jai)] + t1 * alpha;

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
                        var t0 = alpha * A[@intCast(j + jaj)];

                        var i: isize = 0;
                        var ibij: isize = jbj;
                        var icij: isize = jcj;
                        while (i < M) {
                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += t0 * B[@intCast(ibij)];

                            i += 1;
                            ibij += 1;
                            icij += 1;
                        }

                        var k: isize = 0;
                        var iakj: isize = jaj;
                        var jbk: isize = 0;
                        while (k < j) {
                            t0 = alpha * A[@intCast(iakj)];

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * B[@intCast(ibik)];

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
                            t0 = alpha * A[@intCast(iajk)];

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * B[@intCast(ibik)];

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
                        var t0 = alpha * A[@intCast(j + jaj)];

                        var i: isize = 0;
                        var ibij: isize = jbj;
                        var icij: isize = jcj;
                        while (i < M) {
                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += t0 * B[@intCast(ibij)];

                            i += 1;
                            ibij += 1;
                            icij += 1;
                        }

                        var k: isize = 0;
                        var iajk: isize = iaj;
                        var jbk: isize = 0;
                        while (k < j) {
                            t0 = alpha * A[@intCast(iajk)];

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * B[@intCast(ibik)];

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
                            t0 = alpha * A[@intCast(iakj)];

                            i = 0;
                            var ibik: isize = jbk;
                            icij = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * B[@intCast(ibik)];

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
        .cfloat => {
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
                                t1.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                t1.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

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
                            C[@intCast(icij)].re += t0.re * A[@intCast(k)].re - t0.im * A[@intCast(k)].im + t1.re * alpha.re - t1.im * alpha.im;
                            C[@intCast(icij)].im += t0.re * A[@intCast(k)].im + t0.im * A[@intCast(k)].re + t1.re * alpha.im + t1.im * alpha.re;

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
                                t1.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                t1.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

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
                            C[@intCast(icij)].re += t0.re * A[@intCast(k)].re - t0.im * A[@intCast(k)].im + t1.re * alpha.re - t1.im * alpha.im;
                            C[@intCast(icij)].im += t0.re * A[@intCast(k)].im + t0.im * A[@intCast(k)].re + t1.re * alpha.im + t1.im * alpha.re;

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
                        t0.re = alpha.re * A[@intCast(i)].re - alpha.im * A[@intCast(i)].im;
                        t0.im = alpha.im * A[@intCast(i)].re + alpha.re * A[@intCast(i)].im;

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
                            t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                            t0.im = alpha.im * A[@intCast(iajk)].re + alpha.re * A[@intCast(iajk)].im;

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
                        t0.re = alpha.re * A[@intCast(i)].re - alpha.im * A[@intCast(i)].im;
                        t0.im = alpha.im * A[@intCast(i)].re + alpha.re * A[@intCast(i)].im;

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
                            t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                            t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.symm only supports simple types."),
        .unsupported => unreachable,
    }
}

test "symm" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(f64, m * m);
    defer a.free(A);
    const B = try a.alloc(f64, n * n);
    defer a.free(B);
    const C = try a.alloc(f64, m * n);
    defer a.free(C);
    const D = try a.alloc(f64, m * n);
    defer a.free(D);

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
    });
    @memcpy(D.ptr, &[_]f64{
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

    blas.symm(f64, .RowMajor, .Left, .Upper, m, n, alpha, A.ptr, m, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(223, D[0]);
    try std.testing.expectEqual(246, D[1]);
    try std.testing.expectEqual(269, D[2]);
    try std.testing.expectEqual(292, D[3]);
    try std.testing.expectEqual(315, D[4]);
    try std.testing.expectEqual(504, D[5]);
    try std.testing.expectEqual(553, D[6]);
    try std.testing.expectEqual(602, D[7]);
    try std.testing.expectEqual(651, D[8]);
    try std.testing.expectEqual(700, D[9]);
    try std.testing.expectEqual(749, D[10]);
    try std.testing.expectEqual(818, D[11]);
    try std.testing.expectEqual(887, D[12]);
    try std.testing.expectEqual(956, D[13]);
    try std.testing.expectEqual(1025, D[14]);
    try std.testing.expectEqual(928, D[15]);
    try std.testing.expectEqual(1011, D[16]);
    try std.testing.expectEqual(1094, D[17]);
    try std.testing.expectEqual(1177, D[18]);
    try std.testing.expectEqual(1260, D[19]);

    blas.symm(f64, .ColumnMajor, .Left, .Upper, m, n, alpha, A.ptr, m, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(849, D[0]);
    try std.testing.expectEqual(944, D[1]);
    try std.testing.expectEqual(1051, D[2]);
    try std.testing.expectEqual(1176, D[3]);
    try std.testing.expectEqual(1349, D[4]);
    try std.testing.expectEqual(1998, D[5]);
    try std.testing.expectEqual(2263, D[6]);
    try std.testing.expectEqual(2570, D[7]);
    try std.testing.expectEqual(2581, D[8]);
    try std.testing.expectEqual(2866, D[9]);
    try std.testing.expectEqual(3211, D[10]);
    try std.testing.expectEqual(3682, D[11]);
    try std.testing.expectEqual(3513, D[12]);
    try std.testing.expectEqual(3914, D[13]);
    try std.testing.expectEqual(4399, D[14]);
    try std.testing.expectEqual(4476, D[15]);
    try std.testing.expectEqual(4109, D[16]);
    try std.testing.expectEqual(4608, D[17]);
    try std.testing.expectEqual(5215, D[18]);
    try std.testing.expectEqual(5936, D[19]);

    blas.symm(f64, .RowMajor, .Left, .Lower, m, n, alpha, A.ptr, m, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(3223, D[0]);
    try std.testing.expectEqual(3564, D[1]);
    try std.testing.expectEqual(3941, D[2]);
    try std.testing.expectEqual(4372, D[3]);
    try std.testing.expectEqual(4947, D[4]);
    try std.testing.expectEqual(6744, D[5]);
    try std.testing.expectEqual(7609, D[6]);
    try std.testing.expectEqual(8600, D[7]);
    try std.testing.expectEqual(8703, D[8]);
    try std.testing.expectEqual(9628, D[9]);
    try std.testing.expectEqual(10493, D[10]);
    try std.testing.expectEqual(11996, D[11]);
    try std.testing.expectEqual(11579, D[12]);
    try std.testing.expectEqual(12872, D[13]);
    try std.testing.expectEqual(14417, D[14]);
    try std.testing.expectEqual(14464, D[15]);
    try std.testing.expectEqual(13479, D[16]);
    try std.testing.expectEqual(15092, D[17]);
    try std.testing.expectEqual(17029, D[18]);
    try std.testing.expectEqual(19308, D[19]);

    blas.symm(f64, .ColumnMajor, .Left, .Lower, m, n, alpha, A.ptr, m, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(9729, D[0]);
    try std.testing.expectEqual(10826, D[1]);
    try std.testing.expectEqual(12019, D[2]);
    try std.testing.expectEqual(13356, D[3]);
    try std.testing.expectEqual(14981, D[4]);
    try std.testing.expectEqual(20550, D[5]);
    try std.testing.expectEqual(23287, D[6]);
    try std.testing.expectEqual(26360, D[7]);
    try std.testing.expectEqual(26329, D[8]);
    try std.testing.expectEqual(29386, D[9]);
    try std.testing.expectEqual(32203, D[10]);
    try std.testing.expectEqual(36868, D[11]);
    try std.testing.expectEqual(35037, D[12]);
    try std.testing.expectEqual(39302, D[13]);
    try std.testing.expectEqual(44239, D[14]);
    try std.testing.expectEqual(44592, D[15]);
    try std.testing.expectEqual(40817, D[16]);
    try std.testing.expectEqual(46146, D[17]);
    try std.testing.expectEqual(52339, D[18]);
    try std.testing.expectEqual(59444, D[19]);

    blas.symm(f64, .RowMajor, .Right, .Upper, m, n, alpha, B.ptr, n, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(29297, D[0]);
    try std.testing.expectEqual(32730, D[1]);
    try std.testing.expectEqual(36435, D[2]);
    try std.testing.expectEqual(40548, D[3]);
    try std.testing.expectEqual(45493, D[4]);
    try std.testing.expectEqual(61910, D[5]);
    try std.testing.expectEqual(70473, D[6]);
    try std.testing.expectEqual(79988, D[7]);
    try std.testing.expectEqual(80127, D[8]);
    try std.testing.expectEqual(89458, D[9]);
    try std.testing.expectEqual(97019, D[10]);
    try std.testing.expectEqual(111576, D[11]);
    try std.testing.expectEqual(106549, D[12]);
    try std.testing.expectEqual(119706, D[13]);
    try std.testing.expectEqual(134767, D[14]);
    try std.testing.expectEqual(134336, D[15]);
    try std.testing.expectEqual(123783, D[16]);
    try std.testing.expectEqual(140406, D[17]);
    try std.testing.expectEqual(159477, D[18]);
    try std.testing.expectEqual(181132, D[19]);

    blas.symm(f64, .ColumnMajor, .Right, .Upper, m, n, alpha, B.ptr, n, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(89281, D[0]);
    try std.testing.expectEqual(99690, D[1]);
    try std.testing.expectEqual(110915, D[2]);
    try std.testing.expectEqual(123364, D[3]);
    try std.testing.expectEqual(137967, D[4]);
    try std.testing.expectEqual(187346, D[5]);
    try std.testing.expectEqual(213163, D[6]);
    try std.testing.expectEqual(241836, D[7]);
    try std.testing.expectEqual(242007, D[8]);
    try std.testing.expectEqual(270154, D[9]);
    try std.testing.expectEqual(292991, D[10]);
    try std.testing.expectEqual(336816, D[11]);
    try std.testing.expectEqual(321483, D[12]);
    try std.testing.expectEqual(361142, D[13]);
    try std.testing.expectEqual(406513, D[14]);
    try std.testing.expectEqual(405408, D[15]);
    try std.testing.expectEqual(373499, D[16]);
    try std.testing.expectEqual(423598, D[17]);
    try std.testing.expectEqual(481041, D[18]);
    try std.testing.expectEqual(546236, D[19]);

    blas.symm(f64, .RowMajor, .Right, .Lower, m, n, alpha, B.ptr, n, C.ptr, n, beta, D.ptr, n);

    try std.testing.expectEqual(268273, D[0]);
    try std.testing.expectEqual(299538, D[1]);
    try std.testing.expectEqual(333267, D[2]);
    try std.testing.expectEqual(370692, D[3]);
    try std.testing.expectEqual(414611, D[4]);
    try std.testing.expectEqual(563018, D[5]);
    try std.testing.expectEqual(640597, D[6]);
    try std.testing.expectEqual(726800, D[7]);
    try std.testing.expectEqual(727561, D[8]);
    try std.testing.expectEqual(812322, D[9]);
    try std.testing.expectEqual(880503, D[10]);
    try std.testing.expectEqual(1012196, D[11]);
    try std.testing.expectEqual(966511, D[12]);
    try std.testing.expectEqual(1085906, D[13]);
    try std.testing.expectEqual(1222549, D[14]);
    try std.testing.expectEqual(1218304, D[15]);
    try std.testing.expectEqual(1122885, D[16]);
    try std.testing.expectEqual(1273626, D[17]);
    try std.testing.expectEqual(1446543, D[18]);
    try std.testing.expectEqual(1642868, D[19]);

    blas.symm(f64, .ColumnMajor, .Right, .Lower, m, n, alpha, B.ptr, n, C.ptr, m, beta, D.ptr, m);

    try std.testing.expectEqual(805169, D[0]);
    try std.testing.expectEqual(898994, D[1]);
    try std.testing.expectEqual(1000211, D[2]);
    try std.testing.expectEqual(1112516, D[3]);
    try std.testing.expectEqual(1244625, D[4]);
    try std.testing.expectEqual(1689918, D[5]);
    try std.testing.expectEqual(1922727, D[6]);
    try std.testing.expectEqual(2181408, D[7]);
    try std.testing.expectEqual(2183877, D[8]);
    try std.testing.expectEqual(2438266, D[9]);
    try std.testing.expectEqual(2642915, D[10]);
    try std.testing.expectEqual(3038100, D[11]);
    try std.testing.expectEqual(2901057, D[12]);
    try std.testing.expectEqual(3259374, D[13]);
    try std.testing.expectEqual(3669435, D[14]);
    try std.testing.expectEqual(3656832, D[15]);
    try std.testing.expectEqual(3370405, D[16]);
    try std.testing.expectEqual(3822778, D[17]);
    try std.testing.expectEqual(4341679, D[18]);
    try std.testing.expectEqual(4930804, D[19]);

    const gamma = Complex(f64).init(2, 2);
    const delta = Complex(f64).init(3, 3);

    const E = try a.alloc(Complex(f64), m * m);
    defer a.free(E);
    const F = try a.alloc(Complex(f64), n * n);
    defer a.free(F);
    const G = try a.alloc(Complex(f64), m * n);
    defer a.free(G);
    const H = try a.alloc(Complex(f64), m * n);
    defer a.free(H);

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
        Complex(f64).init(16, 16),
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
    @memcpy(G.ptr, &[_]Complex(f64){
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
    @memcpy(H.ptr, &[_]Complex(f64){
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

    blas.symm(Complex(f64), .RowMajor, .Left, .Upper, m, n, gamma, E.ptr, m, G.ptr, n, delta, H.ptr, n);

    try std.testing.expectEqual(-440, H[0].re);
    try std.testing.expectEqual(446, H[0].im);
    try std.testing.expectEqual(-480, H[1].re);
    try std.testing.expectEqual(492, H[1].im);
    try std.testing.expectEqual(-520, H[2].re);
    try std.testing.expectEqual(538, H[2].im);
    try std.testing.expectEqual(-560, H[3].re);
    try std.testing.expectEqual(584, H[3].im);
    try std.testing.expectEqual(-600, H[4].re);
    try std.testing.expectEqual(630, H[4].im);
    try std.testing.expectEqual(-972, H[5].re);
    try std.testing.expectEqual(1008, H[5].im);
    try std.testing.expectEqual(-1064, H[6].re);
    try std.testing.expectEqual(1106, H[6].im);
    try std.testing.expectEqual(-1156, H[7].re);
    try std.testing.expectEqual(1204, H[7].im);
    try std.testing.expectEqual(-1248, H[8].re);
    try std.testing.expectEqual(1302, H[8].im);
    try std.testing.expectEqual(-1340, H[9].re);
    try std.testing.expectEqual(1400, H[9].im);
    try std.testing.expectEqual(-1432, H[10].re);
    try std.testing.expectEqual(1498, H[10].im);
    try std.testing.expectEqual(-1564, H[11].re);
    try std.testing.expectEqual(1636, H[11].im);
    try std.testing.expectEqual(-1696, H[12].re);
    try std.testing.expectEqual(1774, H[12].im);
    try std.testing.expectEqual(-1828, H[13].re);
    try std.testing.expectEqual(1912, H[13].im);
    try std.testing.expectEqual(-1960, H[14].re);
    try std.testing.expectEqual(2050, H[14].im);
    try std.testing.expectEqual(-1760, H[15].re);
    try std.testing.expectEqual(1856, H[15].im);
    try std.testing.expectEqual(-1920, H[16].re);
    try std.testing.expectEqual(2022, H[16].im);
    try std.testing.expectEqual(-2080, H[17].re);
    try std.testing.expectEqual(2188, H[17].im);
    try std.testing.expectEqual(-2240, H[18].re);
    try std.testing.expectEqual(2354, H[18].im);
    try std.testing.expectEqual(-2400, H[19].re);
    try std.testing.expectEqual(2520, H[19].im);

    blas.symm(Complex(f64), .ColumnMajor, .Left, .Upper, m, n, gamma, E.ptr, m, G.ptr, m, delta, H.ptr, m);

    try std.testing.expectEqual(-3018, H[0].re);
    try std.testing.expectEqual(378, H[0].im);
    try std.testing.expectEqual(-3328, H[1].re);
    try std.testing.expectEqual(448, H[1].im);
    try std.testing.expectEqual(-3662, H[2].re);
    try std.testing.expectEqual(542, H[2].im);
    try std.testing.expectEqual(-4032, H[3].re);
    try std.testing.expectEqual(672, H[3].im);
    try std.testing.expectEqual(-4498, H[4].re);
    try std.testing.expectEqual(898, H[4].im);
    try std.testing.expectEqual(-6912, H[5].re);
    try std.testing.expectEqual(1080, H[5].im);
    try std.testing.expectEqual(-7718, H[6].re);
    try std.testing.expectEqual(1334, H[6].im);
    try std.testing.expectEqual(-8608, H[7].re);
    try std.testing.expectEqual(1672, H[7].im);
    try std.testing.expectEqual(-8906, H[8].re);
    try std.testing.expectEqual(1418, H[8].im);
    try std.testing.expectEqual(-9752, H[9].re);
    try std.testing.expectEqual(1712, H[9].im);
    try std.testing.expectEqual(-10718, H[10].re);
    try std.testing.expectEqual(2126, H[10].im);
    try std.testing.expectEqual(-12056, H[11].re);
    try std.testing.expectEqual(2672, H[11].im);
    try std.testing.expectEqual(-12114, H[12].re);
    try std.testing.expectEqual(1938, H[12].im);
    try std.testing.expectEqual(-13312, H[13].re);
    try std.testing.expectEqual(2344, H[13].im);
    try std.testing.expectEqual(-14678, H[14].re);
    try std.testing.expectEqual(2918, H[14].im);
    try std.testing.expectEqual(-14232, H[15].re);
    try std.testing.expectEqual(3672, H[15].im);
    try std.testing.expectEqual(-13978, H[16].re);
    try std.testing.expectEqual(2458, H[16].im);
    try std.testing.expectEqual(-15456, H[17].re);
    try std.testing.expectEqual(2976, H[17].im);
    try std.testing.expectEqual(-17150, H[18].re);
    try std.testing.expectEqual(3710, H[18].im);
    try std.testing.expectEqual(-19072, H[19].re);
    try std.testing.expectEqual(4672, H[19].im);

    blas.symm(Complex(f64), .RowMajor, .Left, .Lower, m, n, gamma, E.ptr, m, G.ptr, n, delta, H.ptr, n);

    try std.testing.expectEqual(-11540, H[0].re);
    try std.testing.expectEqual(-6568, H[0].im);
    try std.testing.expectEqual(-12792, H[1].re);
    try std.testing.expectEqual(-7176, H[1].im);
    try std.testing.expectEqual(-14188, H[2].re);
    try std.testing.expectEqual(-7784, H[2].im);
    try std.testing.expectEqual(-15800, H[3].re);
    try std.testing.expectEqual(-8392, H[3].im);
    try std.testing.expectEqual(-17988, H[4].re);
    try std.testing.expectEqual(-9000, H[4].im);
    try std.testing.expectEqual(-25476, H[5].re);
    try std.testing.expectEqual(-15996, H[5].im);
    try std.testing.expectEqual(-28796, H[6].re);
    try std.testing.expectEqual(-17512, H[6].im);
    try std.testing.expectEqual(-32620, H[7].re);
    try std.testing.expectEqual(-19028, H[7].im);
    try std.testing.expectEqual(-32892, H[8].re);
    try std.testing.expectEqual(-20544, H[8].im);
    try std.testing.expectEqual(-36452, H[9].re);
    try std.testing.expectEqual(-22060, H[9].im);
    try std.testing.expectEqual(-40252, H[10].re);
    try std.testing.expectEqual(-24056, H[10].im);
    try std.testing.expectEqual(-46084, H[11].re);
    try std.testing.expectEqual(-26252, H[11].im);
    try std.testing.expectEqual(-44236, H[12].re);
    try std.testing.expectEqual(-28448, H[12].im);
    try std.testing.expectEqual(-49228, H[13].re);
    try std.testing.expectEqual(-30644, H[13].im);
    try std.testing.expectEqual(-55228, H[14].re);
    try std.testing.expectEqual(-32840, H[14].im);
    try std.testing.expectEqual(-55784, H[15].re);
    try std.testing.expectEqual(-29608, H[15].im);
    try std.testing.expectEqual(-51612, H[16].re);
    try std.testing.expectEqual(-32256, H[16].im);
    try std.testing.expectEqual(-57832, H[17].re);
    try std.testing.expectEqual(-34904, H[17].im);
    try std.testing.expectEqual(-65348, H[18].re);
    try std.testing.expectEqual(-37552, H[18].im);
    try std.testing.expectEqual(-74232, H[19].re);
    try std.testing.expectEqual(-40200, H[19].im);

    blas.symm(Complex(f64), .ColumnMajor, .Left, .Lower, m, n, gamma, E.ptr, m, G.ptr, m, delta, H.ptr, m);

    try std.testing.expectEqual(-15036, H[0].re);
    try std.testing.expectEqual(-54204, H[0].im);
    try std.testing.expectEqual(-17116, H[1].re);
    try std.testing.expectEqual(-59636, H[1].im);
    try std.testing.expectEqual(-19604, H[2].re);
    try std.testing.expectEqual(-65524, H[2].im);
    try std.testing.expectEqual(-22704, H[3].re);
    try std.testing.expectEqual(-72096, H[3].im);
    try std.testing.expectEqual(-27244, H[4].re);
    try std.testing.expectEqual(-80684, H[4].im);
    try std.testing.expectEqual(-29076, H[5].re);
    try std.testing.expectEqual(-123780, H[5].im);
    try std.testing.expectEqual(-34772, H[6].re);
    try std.testing.expectEqual(-138004, H[6].im);
    try std.testing.expectEqual(-41896, H[7].re);
    try std.testing.expectEqual(-153824, H[7].im);
    try std.testing.expectEqual(-37484, H[8].re);
    try std.testing.expectEqual(-159868, H[8].im);
    try std.testing.expectEqual(-44180, H[9].re);
    try std.testing.expectEqual(-174532, H[9].im);
    try std.testing.expectEqual(-50036, H[10].re);
    try std.testing.expectEqual(-191476, H[10].im);
    try std.testing.expectEqual(-61256, H[11].re);
    try std.testing.expectEqual(-215248, H[11].im);
    try std.testing.expectEqual(-47964, H[12].re);
    try std.testing.expectEqual(-217452, H[12].im);
    try std.testing.expectEqual(-57124, H[13].re);
    try std.testing.expectEqual(-238244, H[13].im);
    try std.testing.expectEqual(-69140, H[14].re);
    try std.testing.expectEqual(-262228, H[14].im);
    try std.testing.expectEqual(-80928, H[15].re);
    try std.testing.expectEqual(-253776, H[15].im);
    try std.testing.expectEqual(-58828, H[16].re);
    try std.testing.expectEqual(-250844, H[16].im);
    try std.testing.expectEqual(-70524, H[17].re);
    try std.testing.expectEqual(-276468, H[17].im);
    try std.testing.expectEqual(-85892, H[18].re);
    try std.testing.expectEqual(-306196, H[18].im);
    try std.testing.expectEqual(-105136, H[19].re);
    try std.testing.expectEqual(-340256, H[19].im);

    blas.symm(Complex(f64), .RowMajor, .Right, .Upper, m, n, gamma, F.ptr, n, G.ptr, n, delta, H.ptr, n);

    try std.testing.expectEqual(117284, H[0].re);
    try std.testing.expectEqual(-207500, H[0].im);
    try std.testing.expectEqual(127056, H[1].re);
    try std.testing.expectEqual(-229752, H[1].im);
    try std.testing.expectEqual(137004, H[2].re);
    try std.testing.expectEqual(-254628, H[2].im);
    try std.testing.expectEqual(147216, H[3].re);
    try std.testing.expectEqual(-283440, H[3].im);
    try std.testing.expectEqual(159220, H[4].re);
    try std.testing.expectEqual(-322684, H[4].im);
    try std.testing.expectEqual(283592, H[5].re);
    try std.testing.expectEqual(-458048, H[5].im);
    try std.testing.expectEqual(308472, H[6].re);
    try std.testing.expectEqual(-517104, H[6].im);
    try std.testing.expectEqual(333968, H[7].re);
    try std.testing.expectEqual(-585344, H[7].im);
    try std.testing.expectEqual(364872, H[8].re);
    try std.testing.expectEqual(-589776, H[8].im);
    try std.testing.expectEqual(388456, H[9].re);
    try std.testing.expectEqual(-653536, H[9].im);
    try std.testing.expectEqual(423500, H[10].re);
    try std.testing.expectEqual(-723716, H[10].im);
    try std.testing.expectEqual(460032, H[11].re);
    try std.testing.expectEqual(-827568, H[11].im);
    try std.testing.expectEqual(505588, H[12].re);
    try std.testing.expectEqual(-793372, H[12].im);
    try std.testing.expectEqual(539760, H[13].re);
    try std.testing.expectEqual(-882504, H[13].im);
    try std.testing.expectEqual(575164, H[14].re);
    try std.testing.expectEqual(-990004, H[14].im);
    try std.testing.expectEqual(517424, H[15].re);
    try std.testing.expectEqual(-1002992, H[15].im);
    try std.testing.expectEqual(573384, H[16].re);
    try std.testing.expectEqual(-926352, H[16].im);
    try std.testing.expectEqual(613896, H[17].re);
    try std.testing.expectEqual(-1037040, H[17].im);
    try std.testing.expectEqual(655992, H[18].re);
    try std.testing.expectEqual(-1171344, H[18].im);
    try std.testing.expectEqual(699760, H[19].re);
    try std.testing.expectEqual(-1330576, H[19].im);

    blas.symm(Complex(f64), .ColumnMajor, .Right, .Upper, m, n, gamma, F.ptr, n, G.ptr, m, delta, H.ptr, m);

    try std.testing.expectEqual(971572, H[0].re);
    try std.testing.expectEqual(-267868, H[0].im);
    try std.testing.expectEqual(1067424, H[1].re);
    try std.testing.expectEqual(-305088, H[1].im);
    try std.testing.expectEqual(1171676, H[2].re);
    try std.testing.expectEqual(-349652, H[2].im);
    try std.testing.expectEqual(1288528, H[3].re);
    try std.testing.expectEqual(-405232, H[3].im);
    try std.testing.expectEqual(1442736, H[4].re);
    try std.testing.expectEqual(-487416, H[4].im);
    try std.testing.expectEqual(2221688, H[5].re);
    try std.testing.expectEqual(-520136, H[5].im);
    try std.testing.expectEqual(2473240, H[6].re);
    try std.testing.expectEqual(-622408, H[6].im);
    try std.testing.expectEqual(2754192, H[7].re);
    try std.testing.expectEqual(-750384, H[7].im);
    try std.testing.expectEqual(2860692, H[8].re);
    try std.testing.expectEqual(-671460, H[8].im);
    try std.testing.expectEqual(3122416, H[9].re);
    try std.testing.expectEqual(-791680, H[9].im);
    try std.testing.expectEqual(3437780, H[10].re);
    try std.testing.expectEqual(-896780, H[10].im);
    try std.testing.expectEqual(3858624, H[11].re);
    try std.testing.expectEqual(-1098432, H[11].im);
    try std.testing.expectEqual(3893208, H[12].re);
    try std.testing.expectEqual(-859680, H[12].im);
    try std.testing.expectEqual(4262744, H[13].re);
    try std.testing.expectEqual(-1024184, H[13].im);
    try std.testing.expectEqual(4691080, H[14].re);
    try std.testing.expectEqual(-1240096, H[14].im);
    try std.testing.expectEqual(4556448, H[15].re);
    try std.testing.expectEqual(-1451904, H[15].im);
    try std.testing.expectEqual(4494908, H[16].re);
    try std.testing.expectEqual(-1054604, H[16].im);
    try std.testing.expectEqual(4948048, H[17].re);
    try std.testing.expectEqual(-1264672, H[17].im);
    try std.testing.expectEqual(5476788, H[18].re);
    try std.testing.expectEqual(-1540836, H[18].im);
    try std.testing.expectEqual(6085328, H[19].re);
    try std.testing.expectEqual(-1886768, H[19].im);

    blas.symm(Complex(f64), .RowMajor, .Right, .Lower, m, n, gamma, F.ptr, n, G.ptr, n, delta, H.ptr, n);

    try std.testing.expectEqual(3717460, H[0].re);
    try std.testing.expectEqual(2111972, H[0].im);
    try std.testing.expectEqual(4116600, H[1].re);
    try std.testing.expectEqual(2287944, H[1].im);
    try std.testing.expectEqual(4562940, H[2].re);
    try std.testing.expectEqual(2467116, H[2].im);
    try std.testing.expectEqual(5080080, H[3].re);
    try std.testing.expectEqual(2651088, H[3].im);
    try std.testing.expectEqual(5789036, H[4].re);
    try std.testing.expectEqual(2867380, H[4].im);
    try std.testing.expectEqual(8223512, H[5].re);
    try std.testing.expectEqual(5106616, H[5].im);
    try std.testing.expectEqual(9284728, H[6].re);
    try std.testing.expectEqual(5554712, H[6].im);
    try std.testing.expectEqual(10511144, H[7].re);
    try std.testing.expectEqual(6014008, H[7].im);
    try std.testing.expectEqual(10593376, H[8].re);
    try std.testing.expectEqual(6570776, H[8].im);
    try std.testing.expectEqual(11738568, H[9].re);
    try std.testing.expectEqual(6995928, H[9].im);
    try std.testing.expectEqual(13000620, H[10].re);
    try std.testing.expectEqual(7626060, H[10].im);
    try std.testing.expectEqual(14867672, H[11].re);
    try std.testing.expectEqual(8284072, H[11].im);
    try std.testing.expectEqual(14254540, H[12].re);
    try std.testing.expectEqual(9104708, H[12].im);
    try std.testing.expectEqual(15855824, H[13].re);
    try std.testing.expectEqual(9720640, H[13].im);
    try std.testing.expectEqual(17787508, H[14].re);
    try std.testing.expectEqual(10358972, H[14].im);
    try std.testing.expectEqual(18020896, H[15].re);
    try std.testing.expectEqual(9317792, H[15].im);
    try std.testing.expectEqual(16643760, H[16].re);
    try std.testing.expectEqual(10325688, H[16].im);
    try std.testing.expectEqual(18632496, H[17].re);
    try std.testing.expectEqual(11055792, H[17].im);
    try std.testing.expectEqual(21046032, H[18].re);
    try std.testing.expectEqual(11814696, H[18].im);
    try std.testing.expectEqual(23907968, H[19].re);
    try std.testing.expectEqual(12604000, H[19].im);

    blas.symm(Complex(f64), .ColumnMajor, .Right, .Lower, m, n, gamma, F.ptr, n, G.ptr, m, delta, H.ptr, m);

    try std.testing.expectEqual(4815764, H[0].re);
    try std.testing.expectEqual(17488996, H[0].im);
    try std.testing.expectEqual(5485208, H[1].re);
    try std.testing.expectEqual(19214392, H[1].im);
    try std.testing.expectEqual(6286652, H[2].re);
    try std.testing.expectEqual(21090988, H[2].im);
    try std.testing.expectEqual(7286096, H[3].re);
    try std.testing.expectEqual(23194384, H[3].im);
    try std.testing.expectEqual(8763384, H[4].re);
    try std.testing.expectEqual(25970832, H[4].im);
    try std.testing.expectEqual(9348960, H[5].re);
    try std.testing.expectEqual(39992112, H[5].im);
    try std.testing.expectEqual(11188176, H[6].re);
    try std.testing.expectEqual(44520192, H[6].im);
    try std.testing.expectEqual(13489392, H[7].re);
    try std.testing.expectEqual(49577472, H[7].im);
    try std.testing.expectEqual(12065412, H[8].re);
    try std.testing.expectEqual(51494844, H[8].im);
    try std.testing.expectEqual(14225320, H[9].re);
    try std.testing.expectEqual(56206088, H[9].im);
    try std.testing.expectEqual(16120868, H[10].re);
    try std.testing.expectEqual(61882852, H[10].im);
    try std.testing.expectEqual(19747776, H[11].re);
    try std.testing.expectEqual(69458256, H[11].im);
    try std.testing.expectEqual(15446448, H[12].re);
    try std.testing.expectEqual(70080792, H[12].im);
    try std.testing.expectEqual(18402240, H[13].re);
    try std.testing.expectEqual(76732704, H[13].im);
    try std.testing.expectEqual(22282032, H[14].re);
    try std.testing.expectEqual(84443016, H[14].im);
    try std.testing.expectEqual(26105472, H[15].re);
    try std.testing.expectEqual(82019904, H[15].im);
    try std.testing.expectEqual(18950716, H[16].re);
    try std.testing.expectEqual(80911844, H[16].im);
    try std.testing.expectEqual(22726312, H[17].re);
    try std.testing.expectEqual(89068664, H[17].im);
    try std.testing.expectEqual(27689908, H[18].re);
    try std.testing.expectEqual(98586284, H[18].im);
    try std.testing.expectEqual(33907504, H[19].re);
    try std.testing.expectEqual(109540304, H[19].im);
}
