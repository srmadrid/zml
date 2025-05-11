const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

pub inline fn spr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    const LENX = N;
    const LENY = N;

    switch (numericType) {
        .bool => @compileError("blas.spr2 does not support bool."),
        .int, .float => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)] + t1 * y[@intCast(iy)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    Ap[@intCast(iaij)] += t0 * x[@intCast(jx)] + t1 * y[@intCast(jy)];

                    j += 1;
                    jaj += j;
                    jx += incx;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    Ap[@intCast(jaj)] += t0 * x[@intCast(jx)] + t1 * y[@intCast(jy)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)] + t1 * y[@intCast(iy)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    j += 1;
                    jaj += N - j + 1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .cfloat => @compileError("blas.spr2 does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.spr2 only supports simple types."),
        .unsupported => unreachable,
    }
}
