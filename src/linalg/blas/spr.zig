const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Scalar = types.Scalar;

pub inline fn spr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.spr does not support bool."),
        .int, .float => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < j) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    Ap[@intCast(iaij)] += t0 * x[@intCast(jx)];

                    j += 1;
                    jaj += j;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    Ap[@intCast(jaj)] += t0 * x[@intCast(jx)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += N - j + 1;
                    jx += incx;
                }
            }
        },
        .cfloat => @compileError("blas.spr does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.spr only supports simple types."),
        .unsupported => unreachable,
    }
}
