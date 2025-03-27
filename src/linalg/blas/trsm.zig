const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Side = blas.Side;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;
const Diag = blas.Diag;

pub inline fn trsm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or m <= 0) return;

    var M = m;
    var N = n;
    var SIDE = side;
    var UPLO = uplo;
    const TRANSA = transA;
    if (order == .RowMajor) {
        M = n;
        N = m;
        SIDE = if (side == .Left) .Right else .Left;
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (SIDE == .Left and lda < @max(1, M)) return;
    if (SIDE == .Right and lda < @max(1, N)) return;
    if (ldb < @max(1, M)) return;

    switch (numericType) {
        .bool => @compileError("blas.trsm does not support bool."),
        .int, .float => {
            if (alpha == 0) {
                var j: isize = 0;
                var jbj: isize = 0;
                while (j < N) {
                    var i: isize = 0;
                    var ijbj: isize = jbj;
                    while (i < M) {
                        B[@intCast(ijbj)] = 0;

                        i += 1;
                        ijbj += 1;
                    }

                    j += 1;
                    jbj += ldb;
                }

                return;
            }

            if (SIDE == .Left) {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    B[@intCast(ibkj)] /= A[@intCast(k + jak)];

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0 / A[@intCast(i + jai)];

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    B[@intCast(ibkj)] /= A[@intCast(k + jak)];

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0 / A[@intCast(i + jai)];

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                }
            } else {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] /= A[@intCast(j + jaj)];

                                    i += 1;
                                    ibij += 1;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] /= A[@intCast(k + jak)];

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] /= A[@intCast(j + jaj)];

                                    i += 1;
                                    ibij += 1;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] /= A[@intCast(k + jak)];

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    }
                }
            }
        },
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) {
                var j: isize = 0;
                var jbj: isize = 0;
                while (j < N) {
                    var i: isize = 0;
                    var ijbj: isize = jbj;
                    while (i < M) {
                        B[@intCast(ijbj)].re = 0;
                        B[@intCast(ijbj)].im = 0;

                        i += 1;
                        ijbj += 1;
                    }

                    j += 1;
                    jbj += ldb;
                }

                return;
            }

            if (SIDE == .Left) {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re + temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = -A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re - temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = -A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re + temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = -A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re - temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = -A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                }
            } else {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re + temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = -A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re - temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = -A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re + temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = -A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re - temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = -A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.trsm only supports simple types."),
        .unsupported => unreachable,
    }
}

test trsm {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const alpha = 2;

    const A = try a.alloc(f64, m * m);
    defer a.free(A);
    const B = try a.alloc(f64, n * n);
    defer a.free(B);
    const C = try a.alloc(f64, m * n);
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

    blas.trsm(f64, .RowMajor, .Left, .Upper, .NoTrans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(-4.545454545454546, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-3.4090909090909096, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-2.272727272727272, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-1.1363636363636356, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.4545454545454544, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-0.3409090909090911, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2272727272727272, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.1136363636363639, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-0.13636363636363635, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.09090909090909091, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.045454545454545456, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(2, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(2.125, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(2.25, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(2.375, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(2.5, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Upper, .NoTrans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(-3.0733471074380176, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.4390495867768598, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-0.21952479338842978, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.14204545454545445, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0.8109504132231403, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.04648760330578503, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-0.023243801652892606, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0284090909090909, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.19800275482093715, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0.056129476584022044, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-0.009814049586776863, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-0.017045454545454544, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.21212121212121202, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.030303030303030238, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-0.34090909090909094, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(0.25, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(0.07954545454545454, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(0.011363636363636362, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(0.005681818181818182, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Upper, .NoTrans, .Unit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(-66.17665289256199, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-20.251033057851238, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-7.719352617079888, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-1.5223829201101895, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-81.10261707988981, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(38.04442148760331, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(12.283057851239668, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(4.640151515151514, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(0.8918732782369129, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(52.38498622589532, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-6.019628099173554, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-1.9431818181818181, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.6969696969696967, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.19696969696969685, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-8.181818181818182, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(0.5, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(0.1590909090909091, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(0.022727272727272724, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(0.011363636363636364, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(0.625, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Upper, .NoTrans, .Unit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(1136.149449035809, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-300.2031680440763, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(30.23278236914591, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-3.044765840220379, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-4713.847796143249, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(1092.548898071625, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-114.63842975206607, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(9.280303030303028, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(1152.915289256198, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-303.3829201101928, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(46.256198347107436, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-3.8863636363636362, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-1228.3333333333335, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(299.24242424242425, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-31.363636363636363, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(1, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-696.4772727272727, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(169.8181818181818, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-18.727272727272727, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(1.25, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Upper, .Trans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(2272.298898071618, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-600.4063360881526, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(60.46556473829182, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-6.089531680440758, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-9427.695592286498, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-393.2499999999977, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(161.92263544536215, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-17.06175390266293, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(386.33494031221295, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(3041.4375573921016, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-361.0576634109687, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(59.99889389765413, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-228.966462142082, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-189.7810126053927, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(630.02696385341, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-100.53147695967915, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-62.91856321479263, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(186.36660510059275, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-51.65023687286084, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(363.8411464855162, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Upper, .Trans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(4544.597796143236, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-3987.300275482081, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-82.49511645379376, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-127.0199881041823, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-18855.391184572996, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(15581.74265381083, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(1291.3581267217623, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(473.1995523415973, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(772.6698806244259, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(369.9209519436789, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-1034.1230701707782, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(25.515129064246196, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-457.932924284164, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(318.35043270167245, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(199.81326538614286, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-106.37749855254894, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-125.83712642958525, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(166.98647372485195, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-58.23918846617944, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(56.20888321252531, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Upper, .Trans, .Unit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(9089.195592286473, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-7974.600550964162, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-164.9902329075875, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-254.0399762083646, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-37710.78236914599, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(12985.094123048715, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(18531.91735537185, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(1276.3795704983697, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(2053.419713665581, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(76161.40664217934, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-120231.49177854198, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-105748.58957658197, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-9355.552143334155, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-12975.11720163063, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-419597.8728570451, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(1302327.610991863, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(1152374.4640270064, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(103049.52303510295, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(140173.73023814402, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(4576838.768390115, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Upper, .Trans, .Unit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(18178.391184572945, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-106841.15702479305, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(904476.0691209588, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-12308192.003819143, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-75421.56473829199, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(403078.01193755737, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-3314922.2020202023, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(45063773.963916026, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(4106.839427331162, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(131788.61614770288, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-1595310.6998800933, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(21819733.78042509, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-18711.10428666831, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(67605.28703008029, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-1346848.6774348784, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(22104155.720812466, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(2304748.928054013, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-11317645.594199859, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(92714063.04998876, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-1253071965.9589553, C[19], 0.0000001);

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

    blas.trsm(f64, .RowMajor, .Left, .Lower, .NoTrans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(2, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(4, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(6, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(8, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(10, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(0.3333333333333333, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-2.333333333333333, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-3.6666666666666665, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-5, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(0.06060606060606072, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.4242424242424247, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.6666666666666669, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-0.9090909090909092, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(0.026515151515151335, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-0.07954545454545454, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-0.18560606060606033, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2916666666666668, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-0.39772727272727265, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Lower, .NoTrans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(4, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(20, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-6.555555555555555, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1.464646464646465, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.9154040404040406, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-7.333333333333333, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0.7777777777777777, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(1.5160697887970618, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(0.2846648301193754, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.8484848484848494, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(0.060606060606060844, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(0.027548209366391265, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(0.1644714187327824, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-0.1590909090909091, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-0.008838383838383746, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-0.004017447199265464, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025109044995408384, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Lower, .NoTrans, .Unit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(8, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(40, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-53.111111111111114, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-2.92929292929293, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-1.830808080808081, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-14.666666666666666, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-198.44444444444446, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(462.1432506887053, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(29.862258953168052, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(16.61111111111111, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(146.78787878787878, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(1624.4995408631773, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-6292.264261937558, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-407.2419651056016, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-223.55303030303028, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-1996.492883379247, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-22109.27591253444, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Lower, .NoTrans, .Unit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(16, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-32, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(176, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-1920, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(80, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-266.22222222222223, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(1617.6969696969697, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-17606.247474747477, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-29.333333333333332, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-338.22222222222223, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(3379.842056932966, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-37675.26905417814, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(33.22222222222222, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(227.1313131313131, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(1559.4132231404965, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-33247.42659550047, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-814.4839302112032, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(1181.8617998163459, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-9822.566574839306, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(67455.28839531683, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Lower, .Trans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(1829.1469859580616, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(1849.7299927651184, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(29629.04436583466, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-3297.7484973703977, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-3042.49046551952, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-861.0608595423932, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(11962.150854272753, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-5887.769996799956, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-4.192065280908562, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1096.3151143841992, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(6281.691725493573, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-6711.216430837298, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-195.4133118373822, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(1715.597723098757, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-11214.530844994371, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-4155.928324437558, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-101.8104912764004, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(147.73272497704323, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-1227.8208218549132, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(8431.911049414604, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Lower, .Trans, .NonUnit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(-916.4292420458551, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-5643.3891974570815, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(5836.791952520447, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-412.2185621712997, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-6514.854855740813, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-2779.841907141685, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(2977.8142457950403, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-735.9712495999945, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(469.5009206979198, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1647.0758857264204, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(2057.2916452039176, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-838.9020538546622, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(139.6216879282283, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(2982.190089398676, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-1472.2881093938552, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-519.4910405546948, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-792.063839692552, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(245.8129686478249, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-1373.0462016210668, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(1053.9888811768255, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Left, .Lower, .Trans, .Unit, m, n, alpha, A.ptr, m, C.ptr, n);

    try std.testing.expectApproxEqRel(774415.4767530857, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(774088.3558279412, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-243844.99819981697, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(1771339.7020084802, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-1293538.167217851, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-188007.0797492387, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-192707.69482769084, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(62596.75121444381, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-432173.36678750784, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(312836.5860965208, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(19699.31450704868, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(22084.111083067237, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-7095.14568357829, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(47155.766227429354, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-34564.24265409248, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-1038.9820811093896, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-1584.127679385104, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(491.6259372956498, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-2746.0924032421335, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(2107.977762353651, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Left, .Lower, .Trans, .Unit, m, n, alpha, A.ptr, m, C.ptr, m);

    try std.testing.expectApproxEqRel(-432033640.9122368, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(274205641.3917423, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-42999842.84460316, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(3542679.4040169604, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-21097837.59311446, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(11836599.752684655, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1887737.4188020332, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(125193.50242888762, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-6982494.818742164, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(3706667.65372058, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-490620.0369795163, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(44168.222166134474, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-713870.995200655, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(420285.6532835262, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-44192.91536155961, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-2077.964162218779, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-571463.4549441561, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(359841.1653977362, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-56083.6511029719, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(4215.955524707302, C[19], 0.0000001);

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

    blas.trsm(f64, .RowMajor, .Right, .Upper, .NoTrans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(2, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(12, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1.4285714285714284, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.6593406593406594, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.4164256795835744, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-0.29982648930017347, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(22, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-2.8571428571428568, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-1.3186813186813189, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.8328513591671488, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-0.5996529786003469, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(32, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-4.285714285714286, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-1.978021978021978, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-1.249277038750723, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-0.8994794679005205, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Upper, .NoTrans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(4, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-3.4285714285714284, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(3.4285714285714284, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-0.40816326530612235, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.18838304552590268, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.28384570894692374, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-3.2109623170351913, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(3.7613814756671897, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2656683975365293, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.1706539784528277, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.11337343670605977, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-3.2613366846845273, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(3.788660154189362, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2607484141684015, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-0.11246078447442209, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-0.07034622994733825, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-3.2988800996574277, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Upper, .NoTrans, .Unit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(8, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-16, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(104, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-1344, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(25433.14285714286, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(6.857142857142857, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-14.530612244897958, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(95.29670329670328, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-1231.3745987962297, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(23302.639910003243, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(7.5227629513343794, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-15.576862697741818, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(101.70530477102575, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-1314.0003011934339, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(24866.058591154666, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(7.577320308378724, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-15.67613744509425, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(102.452217066669, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-1323.6957756209272, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(25049.409269128242, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Upper, .NoTrans, .Unit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(16, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-32, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(208, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-2688, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(50770.28571428572, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(205.71428571428572, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1277.061224489796, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(16318.593406593407, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-611882.1777690211, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(44488.70839143506, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(13051.780219780221, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-166286.27460451637, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(10150731.753309065, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-806409.8945053607, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-166817.8859574081, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(2758760.0096098236, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-230661585.62871927, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(18326948.36527407, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(3724518.2733103833, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-62638118.151138686, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Upper, .Trans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(-2323.0588118648275, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-241.58823432206924, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(254.49346443030691, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-4558.339849624061, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(4061.6228571428574, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(39592.01889520284, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(119.59655105236119, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(71801.66819282011, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-68155.06784022834, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(3559.0966713148046, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-930036.3708109302, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-1826895.7113572236, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(1653336.076219747, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-70837.43007783518, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-13345.430876592649, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(138181787.0101009, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-68886331.021212, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(2498756.1679683374, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(5666843.4520232985, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-5011049.452091095, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Upper, .Trans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(1704290.5835872542, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-130213.80358934174, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-43709.94071935251, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-6195211.477108199, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(232686.66604000132, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(2689.906889685504, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(227041.52502223247, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-280637.08774039487, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-139900.64337390647, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(6825.4683770230495, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-150314.50995978684, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-20412801.4154428, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(7135180.405956133, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-259962.45803446724, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-574054.2257704168, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(15051830.998327196, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-5510906.48169696, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(199900.493437467, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(453347.4761618639, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-400883.9561672876, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Upper, .Trans, .Unit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(3564271230.245592, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-2183252630.395963, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(296702434.4787916, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-21697889.595816452, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(465373.33208000264, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(82217244.48393933, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-50494582.75152384, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(6973442.078810807, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-552820.0218287349, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(13650.936754046099, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-3296772222.8035583, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(2030270724.8233414, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-282699430.02144355, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(22442244.114747737, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-1148108.4515408336, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-2644602012.4615192, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(1642616788.2833724, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-224762425.1143197, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(16942053.199015234, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-801767.9123345752, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Upper, .Trans, .Unit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(76126470487030.11, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-10388048004726.357, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(380515810914.96906, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(5563365514742.334, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-14968326407956.166, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(2041847492758.063, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-74659332113.82468, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-1094012962472.1334, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(1353836606656.5269, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-184663557229.77097, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(7306386975.442815, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(98610367745.9717, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-79411004697.64476, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(10833480893.716843, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-815514770.4558129, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-5250719165.130979, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(3285233576.566745, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-449524850.2286394, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(33884106.39803047, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-1603535.8246691504, C[19], 0.0000001);

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

    blas.trsm(f64, .RowMajor, .Right, .Lower, .NoTrans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(-2.0728744939271255, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2591093117408906, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044535, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.08421052631578954, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0.4, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-1.5546558704453437, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-0.19433198380566793, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.09716599190283397, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06315789473684225, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0.8, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-1.0364372469635634, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044543, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522273, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789458, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(1.2, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-0.5182186234817799, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522249, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-0.03238866396761125, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-0.021052631578947666, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(1.6, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Lower, .NoTrans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(-4.359647171494604, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(0.28678485375693974, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-0.1164013740361492, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.14264799338739273, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0.12603480048611096, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.5844683337118878, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(0.11985454136743307, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.008648466385521568, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0022685177596748173, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0.12790244062351455, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2954501794817159, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(0.036222852366044174, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-0.001363733219688919, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0017046665246111133, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(0.1280886426592798, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-0.18928617089281893, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0051821862348178, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025910931174089004, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0016842105263158134, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(0.128, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Lower, .NoTrans, .Unit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(6058.950402352815, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(-1193.2654342801884, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(107.99899181151261, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-6.334966410108112, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(0.2520696009722219, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(5883.418225039397, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-1157.2052707914288, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(104.68856613661215, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-6.143854185448048, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(0.2558048812470291, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(5893.311148191252, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-1159.0069661853172, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(104.83515022373751, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-6.151664180694653, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(0.2561772853185596, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(5889.752342637971, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-1158.2504939271255, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(104.75944939271255, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-6.147368421052632, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(0.256, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Lower, .NoTrans, .Unit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(-7389740.5995342145, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(642423.5044605597, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(118887.57019399949, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(-2007894.2324445434, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(4538871.773068049, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-398066.3518808232, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-72579.87954976199, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(1231044.0711746605, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-616820.4336960121, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(55695.2663849596, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(8521.344069235689, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-167095.39952623384, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(46539.6900575325, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-4202.681304069892, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(246.4070914127424, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(11769.264685275943, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-2316.500987854251, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(209.5188987854251, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-12.294736842105264, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(0.512, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Lower, .Trans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(-14779481.199068429, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(12851676.31476167, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(660919.4271538975, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(109519.14107827547, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(755214.5436559251, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-796132.7037616464, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(661662.3519243364, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(252276.89697195345, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(-225513.77554813487, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(75342.70208917606, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(17042.688138471378, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-62349.56112618513, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(50292.65723893088, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-6653.350735752127, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(689.44036854561, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(23538.529370551885, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(-20837.739742717113, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-650.1467973765024, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-562.9916538708635, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-296.50569637011273, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Lower, .Trans, .NonUnit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(-29558962.398136858, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(25703352.62952334, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(1321838.854307795, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(219038.28215655094, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(8661193.411940794, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-7571281.52379571, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(-188621.8578238453, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(9496.747090115003, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(1456639.4113683042, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1260701.5610789226, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-186342.79415812923, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(-65983.68811099563, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(1052249.602126415, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(-896597.9815946613, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-51556.878835035335, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(485.5178872528302, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(729864.7671495227, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-638510.6062070917, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-35912.887501424746, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-8428.277166270682, C[19], 0.0000001);

    blas.trsm(f64, .RowMajor, .Right, .Lower, .Trans, .Unit, m, n, alpha, B.ptr, n, C.ptr, n);

    try std.testing.expectApproxEqRel(-59117924.796273716, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(406114254.036689, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-4220430197.9726415, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(70010126118.18852, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-1590848847064.4153, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-15142563.04759142, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(90478134.56990083, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-919150427.8211241, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(15251773700.676119, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-346577155516.0015, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-372685.58831625845, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(2104146.1536755594, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-19045713.16837504, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(311222125.86813706, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(-7069847549.746554, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(971.0357745056604, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(1453903.3196520114, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-18734542.441757884, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(312418045.1701628, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(-7099161729.261417, C[19], 0.0000001);

    blas.trsm(f64, .ColumnMajor, .Right, .Lower, .Trans, .Unit, m, n, alpha, B.ptr, n, C.ptr, m);

    try std.testing.expectApproxEqRel(-118235849.59254743, C[0], 0.0000001);
    try std.testing.expectApproxEqRel(812228508.073378, C[1], 0.0000001);
    try std.testing.expectApproxEqRel(-8440860395.945283, C[2], 0.0000001);
    try std.testing.expectApproxEqRel(140020252236.37704, C[3], 0.0000001);
    try std.testing.expectApproxEqRel(-3181461222429.6455, C[4], 0.0000001);
    try std.testing.expectApproxEqRel(-1654742142.2419388, C[5], 0.0000001);
    try std.testing.expectApproxEqRel(17062677061.030367, C[6], 0.0000001);
    try std.testing.expectApproxEqRel(-281878805328.39636, C[7], 0.0000001);
    try std.testing.expectApproxEqRel(25482548034387.293, C[8], 0.0000001);
    try std.testing.expectApproxEqRel(-682353059418.2877, C[9], 0.0000001);
    try std.testing.expectApproxEqRel(-111179580671.58372, C[10], 0.0000001);
    try std.testing.expectApproxEqRel(1834973894210.3472, C[11], 0.0000001);
    try std.testing.expectApproxEqRel(-328122086627583.25, C[12], 0.0000001);
    try std.testing.expectApproxEqRel(9565209041355.648, C[13], 0.0000001);
    try std.testing.expectApproxEqRel(1422573782337.1868, C[14], 0.0000001);
    try std.testing.expectApproxEqRel(-23712806277992.73, C[15], 0.0000001);
    try std.testing.expectApproxEqRel(6212018718347207, C[16], 0.0000001);
    try std.testing.expectApproxEqRel(-181056436126041.47, C[17], 0.0000001);
    try std.testing.expectApproxEqRel(-26911579569210.215, C[18], 0.0000001);
    try std.testing.expectApproxEqRel(448836005615343, C[19], 0.0000001);

    const beta = Complex(f64).init(2, 2);

    const D = try a.alloc(Complex(f64), m * m);
    defer a.free(D);
    const E = try a.alloc(Complex(f64), n * n);
    defer a.free(E);
    const F = try a.alloc(Complex(f64), m * n);
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
        Complex(f64).init(16, 16),
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
    });

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .NoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-4.545454545454546, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.545454545454546, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4090909090909096, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4090909090909096, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.272727272727272, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.272727272727272, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1363636363636356, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1363636363636356, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4545454545454544, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4545454545454544, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3409090909090911, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3409090909090911, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2272727272727272, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2272727272727272, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1136363636363639, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1136363636363639, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.13636363636363635, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.13636363636363635, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.09090909090909091, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.09090909090909091, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.045454545454545456, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.045454545454545456, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.125, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.125, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.25, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.25, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.375, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.375, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .NoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-3.0733471074380168, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.0733471074380168, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4390495867768598, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.43904958677685973, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21952479338842978, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21952479338842978, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14204545454545445, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14204545454545445, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.8109504132231403, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.8109504132231403, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04648760330578503, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04648760330578503, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.023243801652892606, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.023243801652892606, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0284090909090909, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0284090909090909, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.19800275482093715, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.19800275482093715, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.05612947658402204, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.05612947658402203, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.009814049586776862, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.00981404958677686, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.017045454545454544, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.017045454545454544, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21212121212121196, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2121212121212119, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03030303030303026, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.030303030303030293, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.34090909090909094, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.34090909090909094, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.07954545454545454, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.07954545454545454, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.011363636363636362, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.011363636363636362, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.005681818181818182, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.005681818181818182, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .NoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(236.60950413223142, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(320.1962809917356, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(76.60950413223141, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(104.2706611570248, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(31.712121212121197, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(16.985537190082635, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(6.809917355371895, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(8.924931129476589, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(329.5399449035813, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(421.88567493112953, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-160.27479338842977, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7.911157024793384, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-51.38636363636364, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.161157024793387, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13.212121212121206, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(5.462121212121209, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4.484848484848484, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12534435261708277, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-209.54545454545456, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.23002754820936866, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(12, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12.039256198347108, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(3.8181818181818183, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.8863636363636367, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5454545454545453, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3939393939393931, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.27272727272727276, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3939393939393938, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(15, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16.363636363636363, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3181818181818182, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.04545454545454545, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.022727272727272728, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.25, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .NoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(21211.212121212095, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(52846.26170798899, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8317.013774104682, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2557.088154269977, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(564.9490358126724, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-311.199724517906, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4.230027548209387, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(31.46969696969697, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-62344.90909090906, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(14860.093663911844, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4564.705234159777, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8745.584022038563, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(229.27685950413206, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(685.6322314049582, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-37.34848484848483, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15.499999999999993, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(22324.027548209368, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8641.46005509642, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1076.9283746556473, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(3509.3581267217633, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-185.10330578512398, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-229.1694214876033, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(15.40909090909091, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1363636363636367, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1088.3636363636358, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(10905.575757575758, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1197.2121212121212, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1200.2424242424242, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(122.72727272727272, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.7272727272727266, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-963.454545454545, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(6471.09090909091, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-679.1818181818182, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-749.9090909090909, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(74.95454545454545, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.045454545454545456, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.5, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .ConjNoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-111550.53936054764, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(39401.8756991401, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(6061.482093663921, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16715.757575757576, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1877.142791551882, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1338.173094582185, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5969.393250688707, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14751.70684113866, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27380.75642791552, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-123986.55257116615, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(2866.5981300609387, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1560.8479004925277, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-179.54958677685948, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(80.45592286501376, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(2312.789485766759, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(223.56083562901745, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(2625.8901515151515, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(7694.728764921947, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1170.3456152433425, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-384.99024334251607, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(41.93989481592787, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-33.38241923365891, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(882.4462809917358, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(134.1818181818181, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2085.0922865013777, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-105.26859504132224, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(218.23209366391185, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-227.89600550964187, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.8367768595041302, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(22.65495867768595, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-808.8863636363637, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-120.43181818181813, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(93.73863636363635, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-84.89772727272728, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.005681818181818343, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(9.369318181818182, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3125, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .ConjNoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-107283.43204253266, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-233114.52260379278, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(5865.670835420318, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(1496.8855584077728, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2757.7996832999415, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1358.8089025586444, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(1843.9633551423326, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-746.1741563360883, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(250596.1260417535, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-59640.060746932126, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-494.2084081169823, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(992.4205432423405, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(23.478611006761877, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-426.8708599423993, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27.94510445362721, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(289.09868572084486, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16028.038586781919, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(7150.476311132867, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(119.23072542353158, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-409.5094632146661, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(28.941431596202438, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-142.79154429342262, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16.772727272727266, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(110.30578512396697, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-166.5566877452211, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4533.992899031639, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(82.82858022372486, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(72.49235641539362, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4.161697032306536, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.19475488354620552, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.03125, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03125, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(100.48652720385655, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1773.986527203857, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(31.136062327823694, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(31.25030130853994, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7567794421487604, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.05223398760330576, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.0390625, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0390625, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .ConjNoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(248934.44489725732, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-696609.3217568875, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(2048530.973759635, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-202519.91137963725, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-298763.97126220894, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-153908.65943526168, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(140711.97374350822, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-58442.641422362874, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(619398.1943780401, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(385466.7055440453, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3784.617987412314, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(4974.58751046655, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-564141.233564989, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-576950.7238918105, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(23994.963863845074, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(126512.4829910677, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-49192.74037267238, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-19178.753506135738, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(1172.7578941809784, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-611.0106123215628, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(341.96595177925013, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-226.20022539444037, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5077.510330578494, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(85338.4194214876, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(7240.341430837299, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10901.113636363638, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(104.99786083980297, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(308.13464187327827, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10.587903831705484, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6.058884297520661, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(3748.9461088154267, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3347.000000000001, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.228477961432489, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(124.77272727272727, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.6180268595041323, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.409090909090909, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.15625, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .ConjNoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(291941123.286915, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(721670897.9579302, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(57934389.48381765, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-106371612.67769988, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8732429.048264388, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(2601213.2239468317, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(398309.2303317422, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(164538.66464229068, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-844272281.3919019, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(241574065.7962333, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(111508881.06817617, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(77506221.32355231, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1414078.8511770617, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9872932.894377662, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-205035.03825444527, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(301014.8937098256, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-446025181.29779387, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(39008389.8004742, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(48430774.9319997, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(49480777.86721584, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(305786.9521890567, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5120073.633836486, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-180831.85950413218, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(160521.81818181823, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(31330.75982970198, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11712.813506970533, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(85.24259120126874, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(993.6203773269889, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16.558039068369645, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-33.29357625845229, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(13041.287190082645, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-717.6728650137763, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-20.281336088154205, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(291.70213498622593, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-9.792871900826446, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14.054235537190083, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .Trans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(2, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(6, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(6, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(8, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(8, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(10, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(10, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.3333333333333333, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.3333333333333333, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(1, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6666666666666666, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6666666666666666, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3333333333333333, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3333333333333333, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6060606060606062, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6060606060606062, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.4545454545454546, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.4545454545454546, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3030303030303029, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3030303030303029, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.15151515151515163, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.15151515151515163, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.37878787878787884, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.37878787878787884, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.28409090909090906, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.28409090909090906, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.18939393939393956, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.18939393939393956, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.09469696969696975, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.09469696969696975, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .Trans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(4, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.36363636363636365, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.36363636363636365, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.15909090909090906, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.15909090909090906, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(20, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(20, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16.22222222222222, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16.22222222222222, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4343434343434345, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4343434343434331, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6275252525252547, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6275252525252567, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6666666666666666, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6666666666666666, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5555555555555555, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5555555555555555, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.06978879706152433, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.06978879706152427, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06416437098255283, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06416437098255275, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6060606060606059, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6060606060606059, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.45454545454545436, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.45454545454545436, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08264462809917357, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08264462809917357, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.030130853994490406, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.030130853994490434, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5681818181818181, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5681818181818181, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4103535353535353, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4103535353535353, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0746097337006428, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07460973370064283, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03264175849403122, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03264175849403123, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .Trans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(16, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4545454545454546, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6363636363636362, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(80, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-96.88888888888889, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16.000000000000004, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(10.262626262626265, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.909090909090905, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.39898989898988635, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.2727272727272725, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(3.939393939393939, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(160, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-162.22222222222223, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-854.2222222222222, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(406.50137741046825, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(159.8383838383839, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(63.904958677685954, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(18.792929292929177, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(24.358585858585915, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(34.57575757575757, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18.575757575757574, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2015.5555555555557, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-224.7750229568411, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(14161.572084481175, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(5827.881772268136, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-973.1000918273651, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2604.748393021121, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(87.43434343434538, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-493.56060606060527, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-598.6666666666666, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-211.08631772268134, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(19231.5886134068, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(26581.614152892562, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .Trans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-32, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(336, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2941.090909090909, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3202.909090909091, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8022.0000000000055, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(87678.72727272728, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-160, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(160, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1857.7777777777778, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-129.77777777777777, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-17048.080808080806, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-17291.47474747475, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27323.3030303031, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(490896.31313131313, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10.424242424242422, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(5.333333333333333, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(723.2323232323232, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(21.01010101010099, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-9401.85123966942, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8292.047750229569, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(7212.656565656542, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(255502.75941230488, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11.131313131313476, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(86.30303030303018, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(593.4747474747486, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-343.8585858585835, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12077.985307621657, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7653.368227731895, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(71180.61662075233, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(331477.35215794324, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(3263.296602387512, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7155.696969696972, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-50932.97796143252, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(18649.749311294785, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(601281.1698806247, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(356244.38383838383, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2851540.576331501, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13768690.494375579, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .ConjTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-64, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-64, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(672, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(6405.818181818182, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5882.181818181818, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-175357.45454545456, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16044, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-320, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-320, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(64.5925925925926, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(640.5925925925926, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(5753.158249158249, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5906.693602693601, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-165767.37710437708, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7147.040404040429, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(58450.707070707074, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(5344.525252525252, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(99.66329966329964, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(347.74410774410774, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1483.9952138464537, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2099.622784316999, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-50124.87514261068, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(4886.924395469846, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(103725.59810223445, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(6150.324150596895, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(10691.375573921028, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(1082.4793388429755, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1415.3721234382447, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2330.0163063138243, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-42563.96172642404, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(10167.997869535493, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(35603.53935359101, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-471.934419957146, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(1156.8167661463158, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5935.299701561067, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-37935.06955922864, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(75687.02410468322, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(1720054.9510545372, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-354788.9318655742, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .ConjTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(128, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-128, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-330.6666666666667, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(117.33333333333333, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1265.3663911845729, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1162.7548209366391, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(1004.5523415977968, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-23008.431129476587, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(640, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-640, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-746.8641975308642, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(554.8641975308642, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(1229.2753800632586, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(1065.2431384552597, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-125.5594454647453, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21685.09375318845, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10689.050505050502, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(116901.41414141415, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(8791.627384960722, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-97384.62401795735, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1134.9478422255702, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6845.317878678261, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-682.6995780073103, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9618.976855751784, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12300.648301193782, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(207451.1962044689, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(9889.713804713814, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-169312.2049790837, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1497.1572069303145, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15555.270163864225, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1333.8074495937635, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11143.346996613034, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(943.8688399142993, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(71207.07870718202, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(1191.8758672584486, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-58953.62666726958, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15617.056949198146, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11563.416437309072, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(57179.82255678855, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(219576.2436760699, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Upper, .ConjTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(512, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-896, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-426.6666666666667, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(205.2231404958675, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4856.242424242424, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(48025.96694214876, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-44007.75757575758, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(2560, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3627.456790123457, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(640, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(2973.397816549331, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(3650.37037037037, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(32996.13748597083, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-52923.3449647995, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-263217.3480257117, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(396492.17630854, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(207232.5028058362, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-172065.9932659933, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(35336.72897267186, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-39756.93760376958, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-24525.82275294896, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-26750.16074426546, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-315197.6333537395, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(977784.4151107029, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-586574.5885113769, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4660710.479134783, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-219740.7120373111, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(2634653.246588938, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(94513.2369881988, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(844078.035703896, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(427085.9034033466, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(163700.837226089, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7691578.578442829, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14962556.300702155, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(61877042.06700195, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(43925728.444105476, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-29595315.333176825, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-30654587.402474638, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Upper, .ConjTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(1024, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1024, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11178.666666666666, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2645.3333333333335, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(110505.96143250688, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-75210.40220385675, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-178453.93939393928, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(2674315.2066115704, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(5120, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(5120, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-59734.91358024691, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5974.913580246914, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(563584.3264972961, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-524352.4636261606, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(370178.60208142246, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(15526557.436894193, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1319419.0486685033, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(266549.65656565665, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(6022943.952657892, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7859510.507091113, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(27991677.406410716, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(124541985.83312044, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2248601266.057583, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1274620409.0991983, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2585964.0969288847, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(1325173.5635139267, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(14452224.448321603, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-30050258.437506378, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(161618666.77532986, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(414654414.9833978, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8410832603.195365, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3221469069.7598, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(526770.1323545151, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(1181573.4812588713, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(6000237.37645172, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-48582286.502811745, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(446348025.98687273, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(751537549.6747103, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-17392224869.993084, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3942699770.014935, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .NoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(2, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(6, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(6, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(8, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(8, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(10, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(10, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3333333333333333, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3333333333333333, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.333333333333333, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.333333333333333, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.6666666666666665, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.6666666666666665, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.06060606060606072, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.06060606060606072, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18181818181818182, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4242424242424247, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4242424242424247, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6666666666666669, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6666666666666669, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9090909090909092, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9090909090909092, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.02651515151515134, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.02651515151515134, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07954545454545453, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07954545454545453, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18560606060606033, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18560606060606033, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29166666666666674, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29166666666666674, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3977272727272727, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3977272727272727, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .NoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(4, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(20, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(20, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6.555555555555555, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6.555555555555555, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.464646464646465, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.464646464646465, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9154040404040407, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9154040404040407, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7.333333333333333, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7.333333333333333, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.7777777777777777, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.7777777777777777, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.5160697887970618, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.5160697887970618, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2846648301193753, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2846648301193752, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8484848484848494, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8484848484848494, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.060606060606060844, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.060606060606060844, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.027548209366391265, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.02754820936639125, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.16447141873278243, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.16447141873278243, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.15909090909090906, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.15909090909090906, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.008838383838383757, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.008838383838383757, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.004017447199265455, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.004017447199265455, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025109044995408536, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025109044995408536, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .NoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(16, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(80, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(80, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-106.22222222222223, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5.85858585858586, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.6616161616161627, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29.333333333333332, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(400, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-396.8888888888889, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1718.2222222222222, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(124.28650137741053, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-58.5858585858586, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(59.724517906336104, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-36.61616161616163, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(33.22222222222223, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-293.3333333333333, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(293.57575757575756, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7248.888888888889, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-751.0009182736453, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(25238.519742883378, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(24068.804809458215, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(1692.6354453627184, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(64.30394857667591, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(996.3131313131314, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(102.13636363636371, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(8392.969696969696, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(407.01423324150625, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(87351.8751147842, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(118914.78150826445, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .NoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-32, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(128, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-704, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-896, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3072, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(18176, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-160, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(160, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1012.4444444444445, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-52.44444444444446, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6482.505050505051, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6731.717171717171, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10222.333333333314, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(150883.34343434343, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(58.666666666666664, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-58.666666666666664, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(1359.111111111111, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(6.222222222222172, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13507.239669421488, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12745.204775022956, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2384.646464646459, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(304108.9439853076, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-139.67676767676772, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6.787878787878796, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-908.0404040404039, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(293.41414141414157, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4186.927456382004, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11258.00183654729, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-72370.27077594117, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(289456.6693067034, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(3256.662993572085, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(3513.878787878789, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(2302.7851239669435, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11344.184573002756, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-78785.2295684114, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(60578.138659320466, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(1501087.7115472911, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(676267.4326216716, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .ConjNoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-64, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-64, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(256, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1792, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1408, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-36352, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6144, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-320, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-320, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(70.81481481481481, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(390.8148148148148, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(2243.905723905724, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2374.16835016835, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-51787.78114478114, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2234.1111111111045, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(30312.888888888887, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(5139.555555555555, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(264.59259259259255, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(719.7037037037036, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(2305.296491081615, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2758.7843169991934, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-57332.44956451568, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(1515.308233853693, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(45614.85338230793, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(3157.614325068863, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(2132.025711662077, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(189.48760330578594, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(2068.188886106241, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1153.717449981913, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-38353.2620866899, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6749.88651276887, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(51346.51910983109, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(855.8787113559811, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(2512.40652739516, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(426.43193296603204, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6558.819214876032, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9530.90943526171, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-86443.87467695205, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(188347.8333120287, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .ConjNoTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(128, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-128, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-128, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(42.666666666666664, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(302.54545454545456, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(333.57575757575756, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(573.0909090909086, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4783.515151515152, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(640, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-640, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-343.6049382716049, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(236.93827160493825, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(475.77920620344855, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(431.7494133251709, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-65.76804662789618, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6755.75383889399, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10279.111111111106, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(60625.77777777777, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(3186.4691358024693, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-20120.395061728388, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1277.2380015194922, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3311.270447188308, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-170.79882049478783, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9779.350273753476, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6315.228650137724, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(91229.70676461585, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(2041.9136822773146, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29699.227017651257, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(632.7022794941638, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5605.28667255239, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(927.059425871259, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8548.00593875028, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1711.7574227119621, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(102693.03821966218, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(428.44182991531096, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-33393.543897422336, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(1927.0907572047413, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7949.267800616572, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-24775.078791186806, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13820.021090360951, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .ConjNoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(512, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-341.3333333333333, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-170.66666666666669, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-62.060606060606005, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1272.2424242424242, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(10713.21212121212, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8420.848484848488, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(2560, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3721.0864197530864, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(2346.6666666666665, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(2648.0595857565554, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(961.7239057239058, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(7329.062493623097, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-20314.558922558925, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-153271.59595959593, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(196363.63636363635, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(33813.72839506171, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21067.85185185184, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(18313.2144282798, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-60137.59575553516, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12272.732008287232, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4572.941388170031, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-76126.54290378526, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(434256.4431180492, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-388069.3953678197, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3379460.404448525, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-138022.78752800563, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(561910.6336830191, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(658901.7271812834, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1083227.7165215898, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5.459218615811551, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(107855.73514593286, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5138239.945390546, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7352010.234029923, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(55899610.42143152, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(40212670.31747729, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6591950.079331793, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9774569.394471677, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .ConjNoTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(1024, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1024, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4437.333333333333, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1024, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(29416.72727272727, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21472.969696969696, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-21558.303030302995, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(587954.4242424241, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(5120, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(5120, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-32615.506172839505, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2748.8395061728397, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(220203.0911131517, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-201847.09968370575, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(76970.11111111159, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(4799697.963371085, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-699270.4646464646, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(86184.08080808085, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(1335935.9281705946, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1545417.337822671, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(3462530.639447315, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(17729460.462934714, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-250191095.9972544, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-151327841.28251377, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1020765.9720436689, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(716259.8004285279, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(6591794.361391692, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11009111.144577082, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(30434869.1547211, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(118843036.91667487, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1755626965.8834887, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-923554613.2981781, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-215722.38872909735, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(215700.5518546341, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(4427584.25102768, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-25843346.2400084, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(181284279.64139685, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(402826806.09331894, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6831641606.926966, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2450801604.2047186, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .Trans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-10.181818181818182, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10.181818181818182, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10.022727272727272, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10.022727272727272, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-9.863636363636363, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9.863636363636363, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-9.704545454545457, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9.704545454545457, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-9.545454545454547, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9.545454545454547, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4545454545454544, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4545454545454544, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4318181818181817, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4318181818181817, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4090909090909087, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4090909090909087, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.386363636363636, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.386363636363636, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3636363636363633, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3636363636363633, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7272727272727273, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7272727272727273, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7159090909090909, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7159090909090909, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7045454545454546, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7045454545454546, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6931818181818182, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6931818181818182, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6818181818181819, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6818181818181819, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.125, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.125, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.25, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.25, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.375, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.375, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .Trans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-11.75103305785124, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.75103305785124, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1751033057851237, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1751033057851237, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4700413223140492, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4700413223140492, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.213068181818182, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.213068181818182, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-17.840909090909097, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-17.840909090909097, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.17045454545454544, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.17045454545454544, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06818181818181819, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06818181818181818, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1761363636363636, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1761363636363636, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7212465564738286, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7212465564738286, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29485192837465557, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29485192837465557, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03460743801652892, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03460743801652892, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08948863636363637, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08948863636363637, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.015840220385675, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.015840220385675, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1015840220385676, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1015840220385676, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.396694214876033, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.396694214876033, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2.272727272727273, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.272727272727273, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2272727272727273, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2272727272727273, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.09090909090909091, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.09090909090909091, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .Trans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(1112.1880165289258, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1444.6508264462811, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(10058.482954545454, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(13518.066632231406, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1368.970385674931, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1390.0309917355373, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(403.6425619834712, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(573.9569559228651, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(1517.2496556473825, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(1807.564393939394, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-287.38429752066116, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13.297520661157018, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2603.5795454545455, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-123.96590909090925, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-300.633608815427, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(27.20179063360881, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-108.0633608815427, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7.912534435261701, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-373.3677685950413, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.8116391184573173, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(15, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15.138429752066116, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(136.36363636363637, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-136.7215909090909, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(13.636363636363637, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-17.699724517906336, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(5.454545454545455, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5.860881542699725, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(18.75, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-20.33677685950413, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(9.090909090909092, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.9090909090909092, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.36363636363636365, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.25, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .Trans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(317994.11707989, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(636553.249311295, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-278105.1150137742, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-61320.26997245176, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(27507.81267217631, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13856.840220385682, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-340.62878787878776, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(1955.1990358126725, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-196174.57851239672, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-125487.1253443526, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(88724.32920110192, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-28233.555096418742, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3653.541322314048, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(8975.322314049587, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-655.6707988980717, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-546.8636363636364, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(180762.93250688704, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-124871.16666666664, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5419.806473829204, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(86220.64118457299, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6502.3595041322305, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6545.731404958678, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(546.1704545454545, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7159090909090651, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-331.3966942148762, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(3048.0964187327822, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-850.7988980716252, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-861.8126721763085, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(126.17355371900825, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.173553719008261, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-266.90909090909093, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(1417.8181818181818, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-371.6363636363636, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-418.1818181818182, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(59.27272727272727, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.7272727272727273, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.5, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(2.5, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .ConjTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(-1320955.7397111617, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(488902.6230069292, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(137674.51101928364, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-550173.0426997247, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(27149.530970865697, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(56168.287294431924, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-212133.4327573253, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-301852.82786543114, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(394675.0672322397, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-383331.30742340774, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(7427.645212455133, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(31545.203522831616, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2981.25, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1381.33126721763, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(1102.7854161449202, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-115.31826529760416, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(41362.57237665915, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(60511.67964771683, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-28741.15647174221, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1844.8176287670094, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1190.4738918106686, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1181.9062734785875, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(241.80371900826447, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(144.79958677685954, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-625.4803405960431, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(3.0931630353117527, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(156.81718006511397, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-164.79401452541947, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.003146130728774, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(23.366782494365136, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-177.22727272727272, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-33.363636363636374, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(52.272727272727266, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-46.45454545454545, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.09090909090909083, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(7.409090909090909, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3125, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .ConjTrans, .NonUnit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(-1360646.6357071854, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2720914.6244045775, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(193018.73556944352, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(41739.5837153149, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-51374.16512609187, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(33863.56464342904, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(37731.60348317889, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-26516.679094665662, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(787516.5802034878, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(784676.190276165, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10808.951126908958, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(3099.913666150613, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(235.4259214989867, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-692.4252840197619, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(14.414783162200514, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(137.84817701811502, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-122385.52941542902, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(101723.29740524452, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(365.32878481155177, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9834.742101866397, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(234.63744792022402, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(183.47656410081277, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-18.099948347107443, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(30.225464876033058, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-113.23560454279823, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1355.6081219786142, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(59.886034461823336, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(52.061498721247034, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4.282596817157298, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.21648111467795883, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.03125, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03125, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(36.63197314049589, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-389.26833677685954, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(17.054106404958677, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(17.445893595041323, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3897210743801653, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.026084710743801625, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.0390625, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0390625, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Left, .Lower, .ConjTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, n);

    try std.testing.expectApproxEqRel(3012076.729212043, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8241820.25161805, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(2345469.826600592, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(382737.6483821436, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-537728.1938463334, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-135712.05814294334, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(2595238.338273481, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2022653.942537382, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2052.50174666376, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(3342144.7490128726, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-37204.9775029281, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-22792.887485549203, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-211053.50833276418, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-235864.17641099202, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(24027.527435322423, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(53786.95307223248, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-451353.45245783654, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-42558.90124194614, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(20569.25814604219, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18992.173378696807, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(100.44676763882251, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(838.1030240420736, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2294.5692148760363, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(23380.351239669424, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(1461.4986505741113, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3984.4410687453046, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(99.03233594396252, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(222.32998372151263, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11.341905863670513, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5.788481404958679, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(851.8006198347109, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-705.2727272727273, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.783574380165291, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(69, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.831611570247934, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.7272727272727275, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.15625, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Left, .Lower, .ConjTrans, .Unit, m, n, beta, D.ptr, m, F.ptr, m);

    try std.testing.expectApproxEqRel(2768504759.2345862, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2570461255.376733, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(128322583.92671204, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1477630527.6132858, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-125375472.50853387, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(95740508.73781578, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(9235784.561621726, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(1145168.791472198, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-85374196.34962702, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(45925866.45536286, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(31257458.62589077, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(14762185.725221165, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1103699.9807390203, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3475609.116954671, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-59518.85127382011, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(155628.9610151098, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27931426.03799832, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1077761.003617438, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(7234559.3287622025, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(7858289.317919366, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(108664.00980124317, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1120379.7599207705, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-51349.84090909092, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(42171.56404958678, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(9300.421934407943, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6302.071159567121, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(108.97806862767413, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(762.8021186697933, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-17.10684891742367, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-34.260774537258385, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.25, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2873.097107438017, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-42.8739669421484, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12.782024793388445, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(160.2964876033058, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7.708677685950413, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.117768595041323, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3125, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .NoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(2, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(12, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(12, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4285714285714284, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4285714285714284, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6593406593406594, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.6593406593406594, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4164256795835744, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4164256795835744, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29982648930017347, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29982648930017347, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(22, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(22, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.8571428571428568, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.8571428571428568, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3186813186813189, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3186813186813189, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8328513591671488, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8328513591671488, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5996529786003473, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5996529786003473, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4.285714285714286, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.285714285714286, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.978021978021978, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.978021978021978, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.249277038750723, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.249277038750723, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.899479467900521, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.899479467900521, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .NoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(4, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4285714285714284, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4285714285714284, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(3.4285714285714284, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(3.4285714285714284, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.40816326530612235, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.40816326530612235, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18838304552590268, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18838304552590268, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.28384570894692374, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.28384570894692374, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2109623170351913, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2109623170351913, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(3.7613814756671897, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(3.7613814756671897, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2656683975365293, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2656683975365293, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1706539784528277, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1706539784528277, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11337343670605987, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11337343670606007, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2613366846845273, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2613366846845273, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(3.788660154189361, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(3.788660154189361, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2607484141684015, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2607484141684015, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11246078447442201, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11246078447442182, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07034622994733825, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07034622994733825, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2988800996574272, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.2988800996574272, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .NoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(16, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-32, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-464, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-48, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5312, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(7104, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(41520, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-240733.7142857143, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(13.714285714285714, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(27.428571428571427, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29.061224489795915, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-410.77551020408157, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-28.835164835164846, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(4893.613814756672, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(6113.25080240754, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(29625.51840293886, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-213608.22096492286, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(15.045525902668759, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(30.091051805337518, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-31.153725395483637, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-444.821639898563, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-37.31780490064864, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(5214.032798774085, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(6698.880692142779, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(35272.295249847724, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-231103.8243866997, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(15.154640616757444, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(30.309281233514888, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-31.352274890188493, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-447.8285271393547, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-37.56981573478117, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(5249.286517017728, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(6744.063795761146, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(35508.58388778371, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-232664.56989938117, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .NoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-32, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(128, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-832, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1024, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3584, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(24832, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(564891.4285714286, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-398427.4285714286, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-795.4285714285714, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-740.5714285714286, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1039.0204081632653, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(11132.734693877552, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(169732.11930926217, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-128367.2213500785, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11561561.559689589, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1975554.2707656724, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(485717.76445000916, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-350941.405123968, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(143918.97017268447, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-100678.48037676609, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3264493.598357687, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-730108.9008573843, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(156171919.6215319, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(240837232.66930512, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15063950.176689755, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2378084.6408869065, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3666154.0338823213, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1311889.0174573776, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(41005861.45451542, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(70859704.0298464, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2231254133.2577763, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9220528200.965605, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(285215310.298787, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(415539112.5606552, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(51137376.0892487, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(118319402.9531865, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(769357966.6635193, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2594648261.869862, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .ConjNoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-64, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-64, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(18.285714285714285, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(54.857142857142854, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(161.05494505494505, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-146.98901098901098, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2727.754771544245, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-281.4667437825333, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(33965.250850202436, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(45495.538230190876, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1481.142857142857, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1590.857142857143, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3603.9650145772594, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(157.6676384839649, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(21624.8257111805, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(26382.728379940312, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(193414.32080407598, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1236186.152467014, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-138485.68222564916, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(1012231.8104748257, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(201356.96075353218, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(287837.94034536893, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(151071.98288681486, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1014951.8682008729, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-37191316.31415118, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(24584610.797510363, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(27540501.157415517, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-19280486.129794642, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(286639.7986494551, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(728743.2618301857, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-141719408.0596928, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(82011722.90903085, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2674927888.292942, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(614069260.0996414, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1677334085.1758285, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-352935279.14820385, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13759735.570656925, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-42700082.75155509, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(183352826.8063734, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(45439822.40158513, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .ConjNoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(128, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-128, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-109.71428571428571, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(36.57142857142857, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(293.97802197802196, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(322.1098901098901, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(562.9334875650666, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5455.50954308849, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13108.43949434025, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(9814.071671486408, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(548.5714285714286, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(391.83673469387753, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-297.02905840515155, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1305.7984814019799, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8020.436812181572, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(10854.672668698846, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(202174.27529739318, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(20805.21396540887, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-156141.50930381933, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21698.130075531575, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-44257.33001781802, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(31910.791884052258, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(163073.59306036224, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(17838.346055632068, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2767766.477840527, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3943258.8262155256, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(2177049.9548075655, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(2919174.8573075607, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-34763.6177736004, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(838.4294383181548, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8780601.579112746, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14939849.875319233, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-46643057.41576298, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(217752089.87668318, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(26288113.98094342, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-136969547.9337386, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(3490110.8808363914, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1130063.1360907261, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(4650749.13847106, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(28989101.262512892, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .ConjNoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(512.0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.0, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1316.5714285714284, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(877.7142857142858, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1918.593406593405, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14786.109890109889, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(194083.83111625217, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(206380.12261422782, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7860282.778259934, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21326.873297528116, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(313.46938775510205, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1880.816326530612, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2371.032582577772, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6340.348957165283, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(25357.97621332646, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(32720.961893408105, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-380741.7198297706, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(372331.63862064667, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-895726.8190491127, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15495734.803916274, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-152336.24380374057, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-24693.076267531527, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(644529.1341520045, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(106537.54315957054, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3126460.7615287867, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9501047.38278126, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(169249381.39982158, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(103768013.97055194, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5277632009.510003, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(1409919997.5997527, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(12318496.592412975, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-47440902.90886396, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-458545481.95199037, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(461736863.92439425, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(406351486.9994859, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7404343436.412836, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(98092894827.33493, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(101310945508.37964, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3883011905965.481, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(43962663469.83505, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .ConjNoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(1024, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1024, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4388.571428571428, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-877.7142857142853, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(33409.40659340659, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-25735.032967032967, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-24592.582995951292, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(800927.90746096, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15690199.809924811, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15763219.303114925, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(28463.020408163262, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16676.57142857143, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-38107.609009066706, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(337443.87428315124, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4672737.9181502145, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4836965.066527998, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(375912354.639576, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(859413.7558631285, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(29116507.724836364, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-32279867.572461385, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3931739.628251543, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4210087.604484324, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(106652730.42900833, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5607866.2598694675, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6234467300.953724, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(6726939251.002443, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(187787103.0545652, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(1651880745.4390454, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13233762630.517948, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7735851824.431312, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1550045170.4938564, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1940029461.4420104, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-21633709953.530266, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-302439550413.82764, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-28458029950.404182, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-47721195341.19597, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(497015161023.3636, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(266857203229.40155, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7865439883305.429, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7759292005366.944, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .Trans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.4, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.4, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(5.996529786003471, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(5.996529786003471, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.4997108155002891, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.4997108155002891, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.19433198380566802, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.19433198380566802, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.10526315789473684, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.10526315789473684, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.8, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.8, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(11.993059572006942, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(11.993059572006942, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.9994216310005782, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.9994216310005782, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.38866396761133604, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.38866396761133604, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.21052631578947367, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.21052631578947367, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.2, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.2, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(17.98958935801041, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(17.98958935801041, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.4991324465008673, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.4991324465008673, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5829959514170041, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5829959514170041, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3157894736842105, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3157894736842105, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.6, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.6, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .Trans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-0.7541088785495345, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.7541088785495345, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10.392524354887755, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10.392524354887755, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.2382715698352973, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.238271569835299, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3167942501624879, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.3167942501624879, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.07942060070762395, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.07942060070762395, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.499199348640766, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.499199348640766, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.028193431481003, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.028193431481003, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3172886971176215, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3172886971176215, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04288079276359695, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04288079276359695, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.09144880263567672, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.09144880263567672, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.6696759026185835, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.6696759026185835, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.470798359960942, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.470798359960942, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11057928221363122, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.11057928221363122, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03675261027061582, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03675261027061582, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.09440443213296398, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.09440443213296398, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.7319567745274114, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.7319567745274114, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.11993059572006938, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.11993059572006938, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.04663967611336033, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.04663967611336033, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.02526315789473684, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.02526315789473684, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.128, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.128, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .Trans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-747.8176216088552, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7011.175320710461, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(2296.602911915639, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(1497.4995264733445, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-246.87738755171955, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(60.02215568730069, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(6.3536480566099165, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.620825057259868, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3176824028304958, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(904.8728199820631, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5545.515421351351, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(1449.6676124800592, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(1634.6887642472288, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-201.75971414053666, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.354758551849653, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(7.3159042108541374, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7.4874273819085255, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.36579521054270686, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(924.1713489578821, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5718.75931972566, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(1497.001671790813, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(1682.9296840981306, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-207.85980822501597, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.048436881677876, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(7.552354570637119, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7.699365011719582, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.37761772853185593, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(1357.7387160208214, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7652.763925359958, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(1971.4461538461537, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(2286.717779063042, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-277.6252631578947, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8.908178137651818, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(10.24, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10.138947368421054, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.512, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .Trans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(203327852.78360853, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-898978264.876339, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(48023246.6819921, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(69788067.03288147, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3179207.5380235966, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1159152.9699249596, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(20425318.55803308, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-68825434.84861836, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(70493162.93461078, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(93636380.21925196, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10717981.54312264, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1076836.2447765719, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(148333.49926930148, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-404144.51897713315, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(4997543.87287007, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(7640312.572272233, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7147991.121185586, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-711196.85138277, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(493755.14829877566, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-439344.6366431183, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(12205.647708558221, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(24683.565111096017, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-552208.4906505316, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-92282.06084172119, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(219117.28404736536, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-189682.64725944414, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-824.681662049861, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(26651.73124227573, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-974.0941828254847, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-982.2847645429363, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(18070.157282761556, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12590.050418678275, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-630.5432504337768, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(8516.327865818392, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-537.4341700404858, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-573.0668825910931, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(40.757894736842104, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.20210526315789323, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.024, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.024, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .ConjTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(1837865136.0926905, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(374622989.665483, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-19743779.65430528, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(13868223.825190134, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7828766.949517252, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2918688.547545212, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(15129951.476212623, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3786232.8199637416, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7490910.417540157, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(5639453.034768863, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(2691012.920559747, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-21598478.03326328, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(1452878.6376216335, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-797022.3913067445, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1256765.6981365168, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(1578351.8886758685, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(37865.3833756079, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-793999.4989289059, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(35147.57093144947, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(39500.41186390205, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-121143.32498534904, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(315986.09491426393, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6834.437159596679, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-196291.14426664248, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(32201.603873428263, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(33805.41505267372, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2888.164110937903, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.779612188365618, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(78.5827811634349, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-77.92753462603878, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(29983.664718193933, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(36556.90609239867, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2533.99011230153, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-85.89276123154804, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(88.188704800931, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-87.30084985493944, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.06495734072022204, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(4.3765362880886425, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08192, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08192, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .ConjTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-739664036.3145018, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(3688566247.5753145, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-64758071.77171356, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-44104455.33680094, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(4505773.250131819, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18134874.176475734, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(10258332.747759746, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(32413722.14434515, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1820483.5150586092, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2150424.630872308, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(6181411.104430954, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(759607.564023214, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(311057.23876211286, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(447057.61537700024, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-502515.3419255616, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-357446.8407883708, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(127080.73499631655, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1136.1960521238984, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6077.8244421336085, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(5828.122785048275, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-48624.59614321479, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18648.888079142038, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(35526.7695088695, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5421.540522881197, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3567.144431942737, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(3645.7088822355504, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8.318863544460626, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-312.92897553124175, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(8.645158890800408, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(8.265307591194052, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3848.103656378807, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(3156.1835117256774, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(6.871420898523833, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-202.71920898412242, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(6.9840679883951555, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(7.0550963840744805, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3501229030470914, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.005196587257617803, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.006553600000000001, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.006553600000000001, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Upper, .ConjTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-46591872861.871735, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(149397883807.02164, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(52405627586.586586, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-23811665571.687786, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4867373451.56044, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2082935213.7950199, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(101327902.41151792, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(257378080.2539944, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(659882.2316273972, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7941816.291861834, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(345596365.6125879, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(167084787.95871577, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-35821720.48066001, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-138370130.1331866, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7321654.5514756655, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(11201422.182994008, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(738115.033259074, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-209815.96070698104, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-23811.894454363766, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-499.4033141706677, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(88275.58374677284, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-263172.0136793377, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-75741.86954408119, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-27517.29985054629, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(4461.1626565189945, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(16759.799619304402, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-82.3924872904704, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1303.720285446929, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.7597025992127122, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(33.82093296398892, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13729.736281126216, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-316.40229870207645, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(170.9468871726587, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-356.48858994634503, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(19.071609779278468, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(13.50047185851267, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.2349269806094185, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.16556463157894719, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.026214400000000002, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Upper, .ConjTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-392130834808.8844, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(205655479840.47678, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(144143934352.04324, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(53175616627.74654, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4763445071.474978, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10669061863.149858, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-175489293.66420484, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(337438868.5883758, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(14067037.173605854, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(10671470.790542293, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(356501674.8398286, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1025194211.3893824, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(200892560.51948217, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-335773454.0381567, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-44612715.980947316, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(20692017.013987504, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(1238315.1006499277, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1028718.9816723056, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(13837.245774937132, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(34376.29527704355, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(700534.0600414323, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-350364.876953307, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(892134.7596124178, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-183782.67176590548, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-41008.17509414605, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(76664.82918649592, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(611.7810575021845, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4068.270843892015, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(52.43052940895177, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(85.05547575798221, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-26829.18454724828, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-28092.277159656584, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(1054.8709542380075, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-371.08340554737265, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(11.142275841531596, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(65.14416327558227, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.1387246980609427, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.8009832243767314, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.052428800000000005, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.052428800000000005, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .NoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-2.0728744939271255, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.0728744939271255, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25910931174089064, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.25910931174089064, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044535, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044535, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08421052631578954, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08421052631578954, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.4, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.4, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.5546558704453437, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.5546558704453437, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.19433198380566793, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.19433198380566793, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.09716599190283397, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.09716599190283397, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06315789473684225, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06315789473684225, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.8, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.8, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.0364372469635634, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.0364372469635634, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044543, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044543, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522273, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522273, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789458, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789458, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.2, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.2, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5182186234817799, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5182186234817799, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522249, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522249, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03238866396761125, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03238866396761125, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.021052631578947666, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.021052631578947666, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(1.6, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.6, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .NoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-4.359647171494604, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4.359647171494604, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2867848537569395, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2867848537569395, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1164013740361492, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1164013740361492, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14264799338739273, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14264799338739279, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.12603480048611096, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.12603480048611096, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5844683337118878, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5844683337118878, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.11985454136743307, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.11985454136743307, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.008648466385521568, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0086484663855216, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0022685177596748173, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0022685177596748173, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.12790244062351455, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.12790244062351455, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2954501794817159, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.2954501794817159, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.03622285236604418, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.036222852366044195, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.001363733219688919, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.001363733219688919, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0017046665246111133, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0017046665246111133, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.1280886426592798, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.1280886426592798, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18928617089281893, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18928617089281893, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0051821862348178, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0051821862348178, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025910931174089004, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0025910931174089004, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0016842105263158134, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0016842105263158134, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.128, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.128, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .NoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(7087.594060080124, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-55293.301016103025, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(4779.55228301211, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(5232.258652236322, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-434.2517243591696, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7901516169745246, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(12.099340846666653, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12.669932820216225, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5041392019444438, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(8274.776862055485, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54629.76055898075, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(4619.093890315948, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(5294.163406523159, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-430.42714353619954, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.63828512420885, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(12.278634299857398, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12.287708370896096, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5116097624940582, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(8285.485244472131, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54700.59758002918, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(4625.227223524398, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(5301.146655739318, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-431.01292989558925, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.666874067760494, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(12.296509695290862, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12.303328361389307, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5123545706371192, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(8279.566510121454, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54661.072869379925, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(4621.935417004049, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(5297.370170040485, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-430.71326315789474, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-11.665101214574882, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(12.288, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-12.294736842105264, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.512, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .NoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(71263327.97972684, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-295176039.1423535, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(9897085.509463476, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(15423477.260957565, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1615804.1194087057, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(3767272.1246219887, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(24540313.17031194, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-83087430.28939013, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(71331317.10482362, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(93760731.85124782, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7062608.661295219, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-776799.7918643384, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1532263.0969396087, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-436852.4680770754, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(18979321.42091695, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(27759011.285642143, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10792384.259549905, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1009666.3424950411, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(494073.8524165287, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-444178.72691045585, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(125270.71855758988, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-66034.8394079562, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3062914.6395497057, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-444376.4168082742, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(422950.9214915832, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-370640.19296825066, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-884.0084210526286, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(34457.047415299385, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-984.6036565096953, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-982.0152908587257, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(125922.23875900276, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-92763.01271851694, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1350.8695060728733, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(19838.61117408907, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-838.0963238866398, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-884.7567287449392, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(49.165473684210525, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.013473684210527637, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.024, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.024, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .ConjNoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(615791323.3216885, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(122398933.40730885, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3590014.9697005916, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(2028626.4058614585, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12537649.9530898, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3443656.860063276, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(18220813.985956635, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4625026.4474019725, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7500858.548099826, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(5706505.36838589, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(3796924.5776710585, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13097563.776303284, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(7438594.67969187, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5382505.210413145, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4418494.286133551, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(4492074.565982486, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(61395.23838536878, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1185967.911354713, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(35534.29815283647, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(39525.90819332229, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-668879.2897066427, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(5454840.53917384, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(29409.539096754837, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-986657.6623904161, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(62042.057626024136, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(65199.81036963831, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3726.2928520551327, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(6.4432725470190535, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(78.56122326869806, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-78.76829252077562, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(219437.87305926852, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(254250.10398118867, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-5901.513816451437, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-165.20736416687828, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(136.11611433067253, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-136.10205570981333, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.1048961772853181, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(5.278790914127423, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08192, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08192, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .ConjNoTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-241662815.8579167, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1235858093.7192059, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11537237.817191869, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9353354.175767913, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(4411024.280246433, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29252427.67650992, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(11726174.769145396, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(38948020.10874516, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1838575.5282960022, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2154231.4661117685, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(3749110.1596775684, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1078608.789834319, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(2496951.0400100118, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2242917.4606542033, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1455457.6428641102, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1268867.007157985, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(189847.40582543847, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(2421.87065037217, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6080.398949543378, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(5889.008434087438, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-839215.1574464784, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-102913.4120043013, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(180615.3605244579, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-20351.003753906214, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6877.050132733875, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(7027.712492545819, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-12.139464959459994, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-403.7037624757548, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(8.735928973902901, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(8.260769087038927, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-26763.17573865144, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(23098.730378449316, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(13.216589133350285, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-472.121105316115, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(10.888164456785066, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(10.889289146453804, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.4223032731301939, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.00839169418282544, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.006553600000000001, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.006553600000000001, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .ConjNoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-291830812555.1305, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1008044173619.4126, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(121346789629.83855, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54689922235.83183, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-7507379217.198353, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3201204420.2855406, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(122059560.0372167, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(308154610.5025109, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(631311.8756315326, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7985613.988815541, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1979380902.5858054, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1461829132.2763495, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-20636461.82604239, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-311651235.62661964, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13813520.598193694, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(14518323.196107019, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(958569.3695062969, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-180806.25672077277, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-23938.814767261632, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-382.7810309118795, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2425971.443866462, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1882937.6236925707, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(76711.65860852692, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(32608.27796875866, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(1631.1050102301451, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(27786.928577369705, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-55.520586462088886, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1624.7202872261666, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.9503197737279478, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(33.993396121883656, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-97602.58638231173, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(85.60975989721373, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(390.3546957355744, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-832.4387780412953, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(29.800655086036485, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(20.904394825314302, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4905355346260387, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.19867755789473684, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.026214400000000002, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .ConjNoTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-2599794008222.1753, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1432448261030.6196, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(336239085475.8946, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(121618684260.67853, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8693921328.071346, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18866554760.883316, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-277602770.6484252, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(589607149.5963379, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(14051378.98747282, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(7355777.048160821, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1034803880.799585, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(6882310323.232733, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(659677847.2456385, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-604332963.6313244, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-99292342.67244877, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(44127808.525947675, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(1668853.917954709, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1320124.1387092713, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(10998.701964728163, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(31123.069306790603, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1087923.344786657, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8618262.963713478, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(5553979.385952993, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(213816.36647232715, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-83540.02279312507, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(125431.16941850333, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(754.3469946452369, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5032.833333401655, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(53.156690073771685, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(85.78163642280215, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-195378.48943641788, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-195033.95324482903, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(2445.5869475537393, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-884.1681646114419, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(17.792520521444366, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(101.41009982270157, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.5837159534626037, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.378426185041551, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.052428800000000005, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.052428800000000005, F[19].im, 0.0000001);

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
    });

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .Trans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(2, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1428571428571428, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1428571428571428, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1758241758241759, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.1758241758241759, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07403123192596876, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.07403123192596879, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.041457489878542655, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.041457489878542586, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(12, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(12, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8.285714285714285, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-8.285714285714285, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.2747252747252757, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.2747252747252757, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5367264314632737, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.536726431463274, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3005668016194328, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3005668016194319, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(22, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(22, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15.428571428571427, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15.428571428571427, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.373626373626376, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.3736263736263767, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9994216310005771, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.9994216310005755, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5596761133603231, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5596761133603244, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-22.57142857142857, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-22.57142857142857, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4725274725274744, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3.4725274725274744, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4621168305378833, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.4621168305378842, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8187854251012143, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.8187854251012121, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .Trans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(4, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.2857142857142856, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.2857142857142856, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3516483516483518, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3516483516483518, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14806246385193753, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.14806246385193758, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1547021399652977, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.1547021399652977, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4.081632653061224, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(4.081632653061224, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.2668759811616948, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.2668759811616948, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3219036602495252, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3219036602495252, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29506428793878203, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.29506428793878203, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.0305424592758308, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.0305424592758308, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(4.8607656080183546, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(4.8607656080183546, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.1413635525839108, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.1413635525839108, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3275807082327434, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.3275807082327435, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0612150298333365, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06121502983333628, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.4927154474699447, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.4927154474699447, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(5.129919607043591, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(5.129919607043591, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7047302903787025, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.7047302903787025, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18601490245166266, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.18601490245166283, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.062176290483739544, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06217629048373961, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2.7262464312224255, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2.7262464312224255, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .Trans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(0, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(16, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(96, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-105.14285714285714, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2237.714285714286, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-67.69230769230762, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(35896.96703296703, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(41396.15500289185, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(177801.87391555822, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1802150.0531636786, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(16.326530612244895, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(97.95918367346937, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-107.02668759811615, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2280.238618524332, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-72.06940427993052, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(36523.51053457819, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(42233.287230248054, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(183655.70059858018, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1836211.556057684, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(19.443062432073418, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(116.65837459244051, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-125.22382880277615, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2688.712753989792, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-112.39855906171084, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(42572.74704989995, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(50254.28249747025, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(238698.9730614121, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2163652.984116936, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(20.519678428174362, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(123.11807056904618, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-129.936991730561, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2810.9442848853687, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-144.63346838154663, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(44019.97349282626, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(52987.757658537506, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(271415.67063305795, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2260499.061274755, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .Trans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(-32, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(32, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(402.2857142857143, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18.285714285714278, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4340.043956043956, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4610.813186813187, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10998.375939849633, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(154586.24407171778, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(3960031.8541584737, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3248696.358496241, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-873.7959183673469, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-735.3469387755102, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-131.56671899528993, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(17883.579277864992, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(326752.90159464604, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-291880.3523093448, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-57681053.25462905, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5533170.369768211, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(4039580.3908635485, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3293390.568061065, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(143269.97415771044, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-115124.64291752204, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4451828.406790434, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-709761.129586435, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(665186909.8457694, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(878591509.0579754, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(-102677392.74073769, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10248057.174425744, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(1348232.5123521076, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4367907.343932085, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(47483540.03640949, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(71374090.13253784, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(4978223771.50061, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29934468392.303402, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(1738586095.9533896, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(2247324411.2342963, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(-118037854.336329, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(60032566.28542504, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(533647431.10546046, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2304773591.948089, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .ConjTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(-64, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-64, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(60.08163265306122, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(169.79591836734693, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(708.0497524453568, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-770.2799178843135, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-16942.882707502842, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-526.0128405277135, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(275510.35846991756, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(317920.5177758748, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1470.6938775510205, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1747.5918367346937, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6370.1888315765855, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(1460.345368916797, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(49540.41061111903, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(50400.39761279786, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(539967.129296435, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6119272.529903805, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-296103.9931272776, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(6151482.567391173, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(230249.2858350441, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(286539.9483154209, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(5432.363451800734, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1517556.6433533418, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-135367765.27858138, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(103494812.30617763, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(129123134.24269868, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(-107739346.4857217, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(731377.8908796873, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(9417160.195163583, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-142748180.26507568, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(94967080.072819, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(8675060838.02818, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(1340949294.652044, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8232703761.222378, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1050681324.9783664, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(151396800.37863842, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(-296811744.70062613, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(98983353.27480741, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(34450161.826169156, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .ConjTrans, .NonUnit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(128, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-128, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-339.59183673469386, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(120.16326530612243, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(1540.559835768627, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1416.0995048907137, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(1052.0256810554274, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(-33885.765415005684, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-90871.00507882138, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(78753.81670569073, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(596.338192419825, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(385.865889212828, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-857.4014870529783, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2224.653810419228, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-14700.692369672368, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(23836.050293178487, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(997317.6231106294, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(34637.517303641886, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-946670.5434470386, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-45819.64639724438, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-43910.951095517215, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(36465.19258759769, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(242274.0575862228, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6012.798399754562, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11586039.173168141, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14312038.32367184, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(12038319.573858796, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(13625462.636257613, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-958842.7737341853, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(50873.71430564286, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10168310.464381173, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15025850.670712283, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-98569180.00546385, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(705402239.2641274, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(74991682.04834096, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-669489357.5979668, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(24538395.214322414, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(12049772.584917188, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(5238940.86267133, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(19940199.21056, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .RowMajor, .Right, .Lower, .ConjTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, n);

    try std.testing.expectApproxEqRel(512, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3991.5102040816328, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(2633.142857142857, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(10917.328825021137, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-67950.5180533752, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(1111373.2332004546, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(1249526.6623024172, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-55669951.701753154, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1660946.5643371881, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(420.94460641399394, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(1964.408163265306, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-11577.6119713433, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15424.891936052285, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(220717.681096573, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(47459.97629815563, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2480960.699380284, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(5223257.309690487, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-73240943.88382377, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-182849057.96249375, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-160752.2873662298, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-14891.517015839054, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(1550436.538264368, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-402642.10372940806, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-6389453.065209363, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29963679.763608947, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(634479900.2769507, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(506532209.5797132, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-26571750972.140938, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(3651040506.624539, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(9715080.41266222, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(-50388322.270186916, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1363903387.3940341, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(1574286534.6144218, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-588230026.9192653, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-35786136985.689384, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(651837841839.1884, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(584647927362.2424, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-28842851713220.875, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(2358782158950.716, F[19].im, 0.0000001);

    blas.trsm(Complex(f64), .ColumnMajor, .Right, .Lower, .ConjTrans, .Unit, m, n, beta, E.ptr, n, F.ptr, m);

    try std.testing.expectApproxEqRel(1024, F[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1024, F[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13249.30612244898, F[1].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2716.7346938775518, F[1].im, 0.0000001);
    try std.testing.expectApproxEqRel(157735.69375679267, F[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-114066.37845670813, F[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-276306.8582039252, F[3].re, 0.0000001);
    try std.testing.expectApproxEqRel(4721799.791005744, F[3].im, 0.0000001);
    try std.testing.expectApproxEqRel(-108022106.27483194, F[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-114661796.53218068, F[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(28845.15451895044, F[5].re, 0.0000001);
    try std.testing.expectApproxEqRel(-16294.437317784255, F[5].im, 0.0000001);
    try std.testing.expectApproxEqRel(-79644.0706707511, F[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(489599.1366122104, F[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-8544470.456006803, F[7].re, 0.0000001);
    try std.testing.expectApproxEqRel(-9459857.98362988, F[7].im, 0.0000001);
    try std.testing.expectApproxEqRel(1766056642.4379594, F[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(58602115.27941038, F[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(219163720.5421796, F[9].re, 0.0000001);
    try std.testing.expectApproxEqRel(-511850484.6722269, F[9].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3702370.01413271, F[10].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4089827.0503873276, F[10].im, 0.0000001);
    try std.testing.expectApproxEqRel(134604306.00267556, F[11].re, 0.0000001);
    try std.testing.expectApproxEqRel(-5375630.857574469, F[11].im, 0.0000001);
    try std.testing.expectApproxEqRel(-23493927221.383266, F[12].re, 0.0000001);
    try std.testing.expectApproxEqRel(23891414326.87819, F[12].im, 0.0000001);
    try std.testing.expectApproxEqRel(4353460986.923592, F[13].re, 0.0000001);
    try std.testing.expectApproxEqRel(12516587218.755836, F[13].im, 0.0000001);
    try std.testing.expectApproxEqRel(-60340356471.482346, F[14].re, 0.0000001);
    try std.testing.expectApproxEqRel(-45840032513.10192, F[14].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1544737662.4401941, F[15].re, 0.0000001);
    try std.testing.expectApproxEqRel(1866618693.48022, F[15].im, 0.0000001);
    try std.testing.expectApproxEqRel(-38969174531.60578, F[16].re, 0.0000001);
    try std.testing.expectApproxEqRel(-921607849860.8365, F[16].im, 0.0000001);
    try std.testing.expectApproxEqRel(-262614894411.0656, F[17].re, 0.0000001);
    try std.testing.expectApproxEqRel(-225045646850.5849, F[17].im, 0.0000001);
    try std.testing.expectApproxEqRel(2258100173704.3096, F[18].re, 0.0000001);
    try std.testing.expectApproxEqRel(2182966537669.085, F[18].im, 0.0000001);
    try std.testing.expectApproxEqRel(-62411485979271.42, F[19].re, 0.0000001);
    try std.testing.expectApproxEqRel(-53034282373263.8, F[19].im, 0.0000001);
}
