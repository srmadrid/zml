const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;

pub inline fn geru(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    if (order == .RowMajor) {
        M = n;
        N = m;
    }

    if (lda < @max(1, M)) return;

    const LENX = m;
    const LENY = n;

    switch (numericType) {
        .bool => @compileError("blas.geru does not support bool."),
        .int, .float => @compileError("blas.geru does not support integer or float types."),
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) return;

            if (order == .ColumnMajor) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * y[@intCast(jy)].re - alpha.im * y[@intCast(jy)].im, alpha.im * y[@intCast(jy)].re + alpha.re * y[@intCast(jy)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
                        A[@intCast(iaij)].re += x[@intCast(ix)].re * t0.re - x[@intCast(ix)].im * t0.im;
                        A[@intCast(iaij)].im += x[@intCast(ix)].re * t0.im + x[@intCast(ix)].im * t0.re;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        A[@intCast(iaij)].re += y[@intCast(iy)].re * t0.re - y[@intCast(iy)].im * t0.im;
                        A[@intCast(iaij)].im += y[@intCast(iy)].re * t0.im + y[@intCast(iy)].im * t0.re;

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.geru only supports simple types."),
        .unsupported => unreachable,
    }
}
