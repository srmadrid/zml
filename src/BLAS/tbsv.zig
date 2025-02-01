const std = @import("std");
const core = @import("../core/core.zig");
const Transpose = @import("../ndarray/ndarray.zig").Transpose;
const Order = @import("../ndarray/ndarray.zig").Order;
const Uplo = @import("../ndarray/ndarray.zig").Uplo;
const Diag = @import("../ndarray/ndarray.zig").Diag;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn tbsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

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

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.tbsv does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            x[@intCast(jx)] /= A[@intCast(k + jaj)];
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

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
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

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
                        var kx: isize = 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / A[@intCast(iaij)];

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
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            x[@intCast(jx)] /= A[@intCast(jaj)];
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

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
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / A[@intCast(jaj)];

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .Complex => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var temp: T = undefined;
                            if (@abs(A[@intCast(k + jaj)].im) < @abs(A[@intCast(k + jaj)].re)) {
                                const temp1 = A[@intCast(k + jaj)].im / A[@intCast(k + jaj)].re;
                                const temp2 = A[@intCast(k + jaj)].re + temp1 * A[@intCast(k + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = A[@intCast(k + jaj)].re / A[@intCast(k + jaj)].im;
                                const temp2 = A[@intCast(k + jaj)].im + temp1 * A[@intCast(k + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re - A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im + A[@intCast(iaij)].im * t0.re;

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
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re - A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im + A[@intCast(iaij)].im * t0.re;

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
                            var temp: T = undefined;
                            if (@abs(A[@intCast(k + jaj)].im) < @abs(A[@intCast(k + jaj)].re)) {
                                const temp1 = -A[@intCast(k + jaj)].im / A[@intCast(k + jaj)].re;
                                const temp2 = A[@intCast(k + jaj)].re - temp1 * A[@intCast(k + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(k + jaj)].re / A[@intCast(k + jaj)].im;
                                const temp2 = -A[@intCast(k + jaj)].im + temp1 * A[@intCast(k + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re + A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im - A[@intCast(iaij)].im * t0.re;

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
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re + A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im - A[@intCast(iaij)].im * t0.re;

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
                        var kx: isize = 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(A[@intCast(iaij)].im) < @abs(A[@intCast(iaij)].re)) {
                                const temp1 = A[@intCast(iaij)].im / A[@intCast(iaij)].re;
                                const temp2 = A[@intCast(iaij)].re + temp1 * A[@intCast(iaij)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = A[@intCast(iaij)].re / A[@intCast(iaij)].im;
                                const temp2 = A[@intCast(iaij)].im + temp1 * A[@intCast(iaij)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

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
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            if (j >= k) kx += incx;

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
                        var kx: isize = 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(A[@intCast(k + jaj)].im) < @abs(A[@intCast(k + jaj)].re)) {
                                const temp1 = -A[@intCast(k + jaj)].im / A[@intCast(k + jaj)].re;
                                const temp2 = A[@intCast(k + jaj)].re - temp1 * A[@intCast(k + jaj)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(k + jaj)].re / A[@intCast(k + jaj)].im;
                                const temp2 = -A[@intCast(k + jaj)].im + temp1 * A[@intCast(k + jaj)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

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
                            var t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var temp: T = undefined;
                            if (@abs(A[@intCast(jaj)].im) < @abs(A[@intCast(jaj)].re)) {
                                const temp1 = A[@intCast(jaj)].im / A[@intCast(jaj)].re;
                                const temp2 = A[@intCast(jaj)].re + temp1 * A[@intCast(jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = A[@intCast(jaj)].re / A[@intCast(jaj)].im;
                                const temp2 = A[@intCast(jaj)].im + temp1 * A[@intCast(jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re - A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im + A[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re - A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im + A[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

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
                        while (j < N) {
                            var temp: T = undefined;
                            if (@abs(A[@intCast(jaj)].im) < @abs(A[@intCast(jaj)].re)) {
                                const temp1 = -A[@intCast(jaj)].im / A[@intCast(jaj)].re;
                                const temp2 = A[@intCast(jaj)].re - temp1 * A[@intCast(jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(jaj)].re / A[@intCast(jaj)].im;
                                const temp2 = -A[@intCast(jaj)].im + temp1 * A[@intCast(jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re + A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im - A[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re -= A[@intCast(iaij)].re * t0.re + A[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= A[@intCast(iaij)].re * t0.im - A[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

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
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(A[@intCast(jaj)].im) < @abs(A[@intCast(jaj)].re)) {
                                const temp1 = A[@intCast(jaj)].im / A[@intCast(jaj)].re;
                                const temp2 = A[@intCast(jaj)].re + temp1 * A[@intCast(jaj)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = A[@intCast(jaj)].re / A[@intCast(jaj)].im;
                                const temp2 = A[@intCast(jaj)].im + temp1 * A[@intCast(jaj)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

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
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(A[@intCast(jaj)].im) < @abs(A[@intCast(jaj)].re)) {
                                const temp1 = -A[@intCast(jaj)].im / A[@intCast(jaj)].re;
                                const temp2 = A[@intCast(jaj)].re - temp1 * A[@intCast(jaj)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(jaj)].re / A[@intCast(jaj)].im;
                                const temp2 = -A[@intCast(jaj)].im + temp1 * A[@intCast(jaj)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];
                            const I0: isize = j + 1;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = I0;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.tbsv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "tbsv" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 1;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);

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
    @memcpy(x1.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x1.ptr, 2);

    try std.testing.expectApproxEqRel(1.46984126984127, x1[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x1[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x1[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x1[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x1[8], 0.0000001);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Upper, .NoTrans, .Unit, n, k, A.ptr, k + 1, x2.ptr, -2);

    try std.testing.expectApproxEqRel(1749, x2[8], 0.0000001);
    try std.testing.expectApproxEqRel(-874, x2[6], 0.0000001);
    try std.testing.expectApproxEqRel(219, x2[4], 0.0000001);
    try std.testing.expectApproxEqRel(-36, x2[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x2[0], 0.0000001);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x3.ptr, 2);

    try std.testing.expectApproxEqRel(0.8242187499999999, x3[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2161458333333333, x3[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x3[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x3[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x3[8], 0.0000001);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);

    @memcpy(x4.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, k, A.ptr, k + 1, x4.ptr, -2);

    try std.testing.expectApproxEqRel(4345, x4[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1448, x4[6], 0.0000001);
    try std.testing.expectApproxEqRel(290, x4[4], 0.0000001);
    try std.testing.expectApproxEqRel(-41, x4[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x4[0], 0.0000001);

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Upper, .Trans, .NonUnit, n, k, A.ptr, k + 1, x5.ptr, 2);

    try std.testing.expectApproxEqRel(1, x5[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x5[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.6, x5[4], 0.0000001);
    try std.testing.expectApproxEqRel(0.05714285714285716, x5[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619047, x5[8], 0.0000001);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Upper, .Trans, .Unit, n, k, A.ptr, k + 1, x6.ptr, -2);

    try std.testing.expectApproxEqRel(1, x6[8], 0.0000001);
    try std.testing.expectApproxEqRel(0, x6[6], 0.0000001);
    try std.testing.expectApproxEqRel(3, x6[4], 0.0000001);
    try std.testing.expectApproxEqRel(-14, x6[2], 0.0000001);
    try std.testing.expectApproxEqRel(117, x6[0], 0.0000001);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, k, A.ptr, k + 1, x7.ptr, 2);

    try std.testing.expectApproxEqRel(0.5, x7[0], 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x7[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x7[4], 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x7[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x7[8], 0.0000001);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Upper, .Trans, .Unit, n, k, A.ptr, k + 1, x8.ptr, -2);

    try std.testing.expectApproxEqRel(1, x8[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1, x8[6], 0.0000001);
    try std.testing.expectApproxEqRel(8, x8[4], 0.0000001);
    try std.testing.expectApproxEqRel(-52, x8[2], 0.0000001);
    try std.testing.expectApproxEqRel(473, x8[0], 0.0000001);

    const x9 = try a.alloc(f64, 2 * n);
    defer a.free(x9);

    @memcpy(x9.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x9.ptr, 2);

    try std.testing.expectApproxEqRel(0.5, x9[0], 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x9[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x9[4], 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x9[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x9[8], 0.0000001);

    const x10 = try a.alloc(f64, 2 * n);
    defer a.free(x10);

    @memcpy(x10.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Lower, .NoTrans, .Unit, n, k, A.ptr, k + 1, x10.ptr, -2);

    try std.testing.expectApproxEqRel(1, x10[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1, x10[6], 0.0000001);
    try std.testing.expectApproxEqRel(8, x10[4], 0.0000001);
    try std.testing.expectApproxEqRel(-52, x10[2], 0.0000001);
    try std.testing.expectApproxEqRel(473, x10[0], 0.0000001);

    const x11 = try a.alloc(f64, 2 * n);
    defer a.free(x11);

    @memcpy(x11.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x11.ptr, 2);

    try std.testing.expectApproxEqRel(1, x11[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x11[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.6, x11[4], 0.0000001);
    try std.testing.expectApproxEqRel(0.05714285714285716, x11[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619047, x11[8], 0.0000001);

    const x12 = try a.alloc(f64, 2 * n);
    defer a.free(x12);

    @memcpy(x12.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, k, A.ptr, k + 1, x12.ptr, -2);

    try std.testing.expectApproxEqRel(1, x12[8], 0.0000001);
    try std.testing.expectApproxEqRel(0, x12[6], 0.0000001);
    try std.testing.expectApproxEqRel(3, x12[4], 0.0000001);
    try std.testing.expectApproxEqRel(-14, x12[2], 0.0000001);
    try std.testing.expectApproxEqRel(117, x12[0], 0.0000001);

    const x13 = try a.alloc(f64, 2 * n);
    defer a.free(x13);

    @memcpy(x13.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Lower, .Trans, .NonUnit, n, k, A.ptr, k + 1, x13.ptr, 2);

    try std.testing.expectApproxEqRel(0.8242187499999999, x13[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.2161458333333333, x13[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x13[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x13[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x13[8], 0.0000001);

    const x14 = try a.alloc(f64, 2 * n);
    defer a.free(x14);

    @memcpy(x14.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .RowMajor, .Lower, .Trans, .Unit, n, k, A.ptr, k + 1, x14.ptr, -2);

    try std.testing.expectApproxEqRel(4345, x14[8], 0.0000001);
    try std.testing.expectApproxEqRel(-1448, x14[6], 0.0000001);
    try std.testing.expectApproxEqRel(290, x14[4], 0.0000001);
    try std.testing.expectApproxEqRel(-41, x14[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x14[0], 0.0000001);

    const x15 = try a.alloc(f64, 2 * n);
    defer a.free(x15);

    @memcpy(x15.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, k, A.ptr, k + 1, x15.ptr, 2);

    try std.testing.expectApproxEqRel(1.46984126984127, x15[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x15[2], 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x15[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x15[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x15[8], 0.0000001);

    const x16 = try a.alloc(f64, 2 * n);
    defer a.free(x16);

    @memcpy(x16.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbsv(f64, .ColumnMajor, .Lower, .Trans, .Unit, n, k, A.ptr, k + 1, x16.ptr, -2);

    try std.testing.expectApproxEqRel(1749, x16[8], 0.0000001);
    try std.testing.expectApproxEqRel(-874, x16[6], 0.0000001);
    try std.testing.expectApproxEqRel(219, x16[4], 0.0000001);
    try std.testing.expectApproxEqRel(-36, x16[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x16[0], 0.0000001);

    const B = try a.alloc(Complex(f64), n * n);
    defer a.free(B);
    const x17 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x17);

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
    @memcpy(x17.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x17.ptr, 2);

    try std.testing.expectApproxEqRel(1.46984126984127, x17[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x17[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x17[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x17[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x17[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[8].im, 0.0000001);

    const x18 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x18);

    @memcpy(x18.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .NoTrans, .Unit, n, k, B.ptr, k + 1, x18.ptr, -2);

    try std.testing.expectApproxEqRel(-6959, x18[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7639, x18[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(3650, x18[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(170, x18[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-477, x18[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(435, x18[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x18[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-76, x18[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x18[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x18[0].im, 0.0000001);

    const x19 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x19);

    @memcpy(x19.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x19.ptr, 2);

    try std.testing.expectApproxEqRel(0.8242187499999999, x19[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21614583333333326, x19[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x19[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x19[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x19[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[8].im, 0.0000001);

    const x20 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x20);

    @memcpy(x20.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .NoTrans, .Unit, n, k, B.ptr, k + 1, x20.ptr, -2);

    try std.testing.expectApproxEqRel(-17309, x20[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18821, x20[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(6022, x20[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(252, x20[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-627, x20[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(577, x20[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x20[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-86, x20[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x20[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x20[0].im, 0.0000001);

    const x21 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x21);

    @memcpy(x21.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x21.ptr, 2);

    try std.testing.expectApproxEqRel(0, x21[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.46984126984127, x21[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x21[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x21[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x21[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x21[8].im, 0.0000001);

    const x22 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x22);

    @memcpy(x22.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x22.ptr, -2);

    try std.testing.expectApproxEqRel(-7639, x22[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6959, x22[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(170, x22[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(3650, x22[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(435, x22[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-477, x22[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-76, x22[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x22[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x22[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x22[0].im, 0.0000001);

    const x23 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x23);

    @memcpy(x23.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x23.ptr, 2);

    try std.testing.expectApproxEqRel(0, x23[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.8242187499999999, x23[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21614583333333326, x23[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x23[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x23[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x23[8].im, 0.0000001);

    const x24 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x24);

    @memcpy(x24.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x24.ptr, -2);

    try std.testing.expectApproxEqRel(-18821, x24[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-17309, x24[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(252, x24[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(6022, x24[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(577, x24[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-627, x24[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-86, x24[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x24[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x24[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x24[0].im, 0.0000001);

    const x25 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x25);

    @memcpy(x25.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .Trans, .NonUnit, n, k, B.ptr, k + 1, x25.ptr, 2);

    try std.testing.expectApproxEqRel(1, x25[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6000000000000001, x25[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.057142857142857065, x25[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619049, x25[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[8].im, 0.0000001);

    const x26 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x26);

    @memcpy(x26.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .Trans, .Unit, n, k, B.ptr, k + 1, x26.ptr, -2);

    try std.testing.expectApproxEqRel(1, x26[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x26[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x26[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x26[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13, x26[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(3, x26[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(100, x26[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(64, x26[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-283, x26[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1307, x26[0].im, 0.0000001);

    const x27 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x27);

    @memcpy(x27.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .Trans, .NonUnit, n, k, B.ptr, k + 1, x27.ptr, 2);

    try std.testing.expectApproxEqRel(0.5, x27[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x27[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x27[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x27[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x27[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[8].im, 0.0000001);

    const x28 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x28);

    @memcpy(x28.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .Trans, .Unit, n, k, B.ptr, k + 1, x28.ptr, -2);

    try std.testing.expectApproxEqRel(1, x28[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x28[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x28[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4, x28[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27, x28[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(13, x28[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(284, x28[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(102, x28[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1633, x28[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3469, x28[0].im, 0.0000001);

    const x29 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x29);

    @memcpy(x29.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x29.ptr, 2);

    try std.testing.expectApproxEqRel(0, x29[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x29[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6000000000000001, x29[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.057142857142857065, x29[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619049, x29[8].im, 0.0000001);

    const x30 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x30);

    @memcpy(x30.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Upper, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x30.ptr, -2);

    try std.testing.expectApproxEqRel(1, x30[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x30[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x30[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x30[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(3, x30[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13, x30[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(64, x30[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(100, x30[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1307, x30[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-283, x30[0].im, 0.0000001);

    const x31 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x31);

    @memcpy(x31.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x31.ptr, 2);

    try std.testing.expectApproxEqRel(0, x31[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x31[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x31[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x31[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x31[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x31[8].im, 0.0000001);

    const x32 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x32);

    @memcpy(x32.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x32.ptr, -2);

    try std.testing.expectApproxEqRel(1, x32[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x32[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4, x32[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x32[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(13, x32[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-27, x32[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(102, x32[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(284, x32[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3469, x32[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1633, x32[0].im, 0.0000001);

    const x33 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x33);

    @memcpy(x33.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x33.ptr, 2);

    try std.testing.expectApproxEqRel(0.5, x33[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x33[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x33[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x33[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x33[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[8].im, 0.0000001);

    const x34 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x34);

    @memcpy(x34.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .NoTrans, .Unit, n, k, B.ptr, k + 1, x34.ptr, -2);

    try std.testing.expectApproxEqRel(1, x34[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x34[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x34[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4, x34[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-27, x34[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(13, x34[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(284, x34[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(102, x34[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1633, x34[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3469, x34[0].im, 0.0000001);

    const x35 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x35);

    @memcpy(x35.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x35.ptr, 2);

    try std.testing.expectApproxEqRel(1, x35[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6000000000000001, x35[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.057142857142857065, x35[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619049, x35[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[8].im, 0.0000001);

    const x36 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x36);

    @memcpy(x36.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .NoTrans, .Unit, n, k, B.ptr, k + 1, x36.ptr, -2);

    try std.testing.expectApproxEqRel(1, x36[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x36[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x36[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x36[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-13, x36[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(3, x36[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(100, x36[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(64, x36[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-283, x36[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1307, x36[0].im, 0.0000001);

    const x37 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x37);

    @memcpy(x37.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x37.ptr, 2);

    try std.testing.expectApproxEqRel(0, x37[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x37[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.125, x37[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.3958333333333333, x37[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.15364583333333337, x37[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.36171875, x37[8].im, 0.0000001);

    const x38 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x38);

    @memcpy(x38.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x38.ptr, -2);

    try std.testing.expectApproxEqRel(1, x38[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x38[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4, x38[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x38[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(13, x38[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-27, x38[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(102, x38[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(284, x38[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3469, x38[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1633, x38[0].im, 0.0000001);

    const x39 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x39);

    @memcpy(x39.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x39.ptr, 2);

    try std.testing.expectApproxEqRel(0, x39[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x39[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6000000000000001, x39[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.057142857142857065, x39[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5047619047619049, x39[8].im, 0.0000001);

    const x40 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x40);

    @memcpy(x40.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x40.ptr, -2);

    try std.testing.expectApproxEqRel(1, x40[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x40[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x40[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x40[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(3, x40[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-13, x40[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(64, x40[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(100, x40[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-1307, x40[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-283, x40[0].im, 0.0000001);

    const x41 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x41);

    @memcpy(x41.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .Trans, .NonUnit, n, k, B.ptr, k + 1, x41.ptr, 2);

    try std.testing.expectApproxEqRel(0.8242187499999999, x41[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21614583333333326, x41[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x41[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x41[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x41[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[8].im, 0.0000001);

    const x42 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x42);

    @memcpy(x42.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .Trans, .Unit, n, k, B.ptr, k + 1, x42.ptr, -2);

    try std.testing.expectApproxEqRel(-17309, x42[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-18821, x42[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(6022, x42[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(252, x42[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-627, x42[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(577, x42[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x42[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-86, x42[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x42[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x42[0].im, 0.0000001);

    const x43 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x43);

    @memcpy(x43.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .Trans, .NonUnit, n, k, B.ptr, k + 1, x43.ptr, 2);

    try std.testing.expectApproxEqRel(1.46984126984127, x43[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x43[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x43[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x43[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x43[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[8].im, 0.0000001);

    const x44 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x44);

    @memcpy(x44.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .Trans, .Unit, n, k, B.ptr, k + 1, x44.ptr, -2);

    try std.testing.expectApproxEqRel(-6959, x44[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-7639, x44[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(3650, x44[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(170, x44[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-477, x44[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(435, x44[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x44[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-76, x44[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x44[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x44[0].im, 0.0000001);

    const x45 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x45);

    @memcpy(x45.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x45.ptr, 2);

    try std.testing.expectApproxEqRel(0, x45[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.8242187499999999, x45[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.21614583333333326, x45[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5729166666666666, x45[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.0625, x45[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5, x45[8].im, 0.0000001);

    const x46 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x46);

    @memcpy(x46.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .RowMajor, .Lower, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x46.ptr, -2);

    try std.testing.expectApproxEqRel(-18821, x46[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-17309, x46[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(252, x46[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(6022, x46[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(577, x46[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-627, x46[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-86, x46[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x46[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x46[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x46[0].im, 0.0000001);

    const x47 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x47);

    @memcpy(x47.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x47.ptr, 2);

    try std.testing.expectApproxEqRel(0, x47[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1.46984126984127, x47[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.23492063492063497, x47[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.6761904761904762, x47[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06349206349206352, x47[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.5555555555555556, x47[8].im, 0.0000001);

    const x48 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x48);

    @memcpy(x48.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbsv(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x48.ptr, -2);

    try std.testing.expectApproxEqRel(-7639, x48[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-6959, x48[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(170, x48[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(3650, x48[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(435, x48[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-477, x48[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-76, x48[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x48[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x48[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x48[0].im, 0.0000001);
}
