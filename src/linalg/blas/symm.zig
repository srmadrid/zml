const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Side = blas.Side;
const Uplo = blas.Uplo;

pub inline fn symm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

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
