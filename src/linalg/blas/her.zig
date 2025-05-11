const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Scalar = types.Scalar;

pub inline fn her(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: Scalar(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.her does not support bool."),
        .int, .float => @compileError("blas.her does not support int or float."),
        .cfloat => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha * x[@intCast(jx)].re, alpha * x[@intCast(jx)].im * (-conj));

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < j) {
                        A[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        A[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    A[@intCast(iaij)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    A[@intCast(iaij)].im = 0;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha * x[@intCast(jx)].re, alpha * x[@intCast(jx)].im * (-conj));
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    A[@intCast(jaj)].im = 0;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        A[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        A[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.her only supports simple types."),
        .unsupported => unreachable,
    }
}
