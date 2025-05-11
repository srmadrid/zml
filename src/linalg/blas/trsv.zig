const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Transpose = blas.Transpose;
const Order = blas.Order;
const Uplo = blas.Uplo;
const Diag = blas.Diag;

pub inline fn trsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var TRANSA = transA;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANSA = if (transA == .NoTrans) .Trans else if (transA == .ConjNoTrans) .ConjTrans else if (transA == .Trans) .NoTrans else .ConjNoTrans;
    }

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.trsv does not support bool."),
        .int, .float => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            x[@intCast(jx)] /= A[@intCast(j + jaj)];
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / A[@intCast(iaij)];

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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

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
            } else {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
                        while (j < N) {
                            x[@intCast(jx)] /= A[@intCast(jaj)];
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)] -= t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / A[@intCast(jaj)];

                            j -= 1;
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    } else {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0 -= A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .cfloat => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var temp: T = undefined;
                            if (@abs(A[@intCast(j + jaj)].im) < @abs(A[@intCast(j + jaj)].re)) {
                                const temp1 = A[@intCast(j + jaj)].im / A[@intCast(j + jaj)].re;
                                const temp2 = A[@intCast(j + jaj)].re + temp1 * A[@intCast(j + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = A[@intCast(j + jaj)].re / A[@intCast(j + jaj)].im;
                                const temp2 = A[@intCast(j + jaj)].im + temp1 * A[@intCast(j + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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
                            if (@abs(A[@intCast(j + jaj)].im) < @abs(A[@intCast(j + jaj)].re)) {
                                const temp1 = -A[@intCast(j + jaj)].im / A[@intCast(j + jaj)].re;
                                const temp2 = A[@intCast(j + jaj)].re - temp1 * A[@intCast(j + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(j + jaj)].re / A[@intCast(j + jaj)].im;
                                const temp2 = -A[@intCast(j + jaj)].im + temp1 * A[@intCast(j + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(A[@intCast(iaij)].im) < @abs(A[@intCast(iaij)].re)) {
                                const temp1 = -A[@intCast(iaij)].im / A[@intCast(iaij)].re;
                                const temp2 = A[@intCast(iaij)].re - temp1 * A[@intCast(iaij)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = -A[@intCast(iaij)].re / A[@intCast(iaij)].im;
                                const temp2 = -A[@intCast(iaij)].im + temp1 * A[@intCast(iaij)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

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

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

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
            } else {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
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

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * A[@intCast(iaij)].re + t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * A[@intCast(iaij)].re + t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
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

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        const ldap1: isize = lda + 1;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += ldap1;
                            jx += incx;
                        }
                    }
                } else if (TRANSA == .Trans) {
                    if (diag == .NonUnit) {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
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
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    } else {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
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
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    } else {
                        const ldap1: isize = lda + 1;
                        var j: isize = N - 1;
                        var jaj: isize = ldap1 * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= ldap1;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.trsv only supports simple types."),
        .unsupported => unreachable,
    }
}
