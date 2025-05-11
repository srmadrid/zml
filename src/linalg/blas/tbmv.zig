const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Transpose = blas.Transpose;
const Order = blas.Order;
const Uplo = blas.Uplo;
const Diag = blas.Diag;

pub inline fn tbmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or k < 0) return;

    const N = n;
    var UPLO = uplo;
    var TRANSA = transA;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANSA = if (transA == .NoTrans) .Trans else if (transA == .ConjNoTrans) .ConjTrans else if (transA == .Trans) .NoTrans else .ConjNoTrans;
    }

    if (lda < k + 1) return;

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.tbmv does not support bool."),
        .int, .float => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] *= A[@intCast(iaij)];

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0: T = 0;
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 + A[@intCast(iaij)] * x[@intCast(jx)];

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0: T = 0;
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] += t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)] *= A[@intCast(jaj)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            t0 *= A[@intCast(jaj)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            }
        },
        .cfloat => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re = t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                            x[@intCast(jx)].im = t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re = t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                            x[@intCast(jx)].im = t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else if (TRANSA == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            t0.re += A[@intCast(iaij)].re * x[@intCast(jx)].re - A[@intCast(iaij)].im * x[@intCast(jx)].im;
                            t0.im += A[@intCast(iaij)].re * x[@intCast(jx)].im + A[@intCast(iaij)].im * x[@intCast(jx)].re;
                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re += t0.re;
                            x[@intCast(jx)].im += t0.im;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            t0.re += A[@intCast(iaij)].re * x[@intCast(jx)].re + A[@intCast(iaij)].im * x[@intCast(jx)].im;
                            t0.im += A[@intCast(iaij)].re * x[@intCast(jx)].im - A[@intCast(iaij)].im * x[@intCast(jx)].re;
                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re += t0.re;
                            x[@intCast(jx)].im += t0.im;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)].re = t0.re * A[@intCast(jaj)].re - t0.im * A[@intCast(jaj)].im;
                            x[@intCast(jx)].im = t0.re * A[@intCast(jaj)].im + t0.im * A[@intCast(jaj)].re;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)].re = t0.re * A[@intCast(jaj)].re + t0.im * A[@intCast(jaj)].im;
                            x[@intCast(jx)].im = t0.im * A[@intCast(jaj)].re - t0.re * A[@intCast(jaj)].im;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANSA == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const temp = t0;
                            t0.re = temp.re * A[@intCast(jaj)].re - temp.im * A[@intCast(jaj)].im;
                            t0.im = temp.re * A[@intCast(jaj)].im + temp.im * A[@intCast(jaj)].re;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const temp = t0;
                            t0.re = temp.re * A[@intCast(jaj)].re + temp.im * A[@intCast(jaj)].im;
                            t0.im = temp.im * A[@intCast(jaj)].re - temp.re * A[@intCast(jaj)].im;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.tbmv only supports simple types."),
        .unsupported => unreachable,
    }
}
