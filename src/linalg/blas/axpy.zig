const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn axpy(comptime T: type, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    switch (numericType) {
        .bool => @compileError("blas.axpy does not support bool."),
        .int, .float => {
            if (alpha == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                const incy2 = incy * 2;
                const incy3 = incy * 3;
                const incy4 = incy * 4;
                while (ix != StX) {
                    y[@intCast(iy)] += alpha * x[@intCast(ix)];
                    y[@intCast(iy + incy)] += alpha * x[@intCast(ix + incx)];
                    y[@intCast(iy + incy2)] += alpha * x[@intCast(ix + incx2)];
                    y[@intCast(iy + incy3)] += alpha * x[@intCast(ix + incx3)];

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)] += alpha * x[@intCast(ix)];

                ix += incx;
                iy += incy;
            }
        },
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    y[@intCast(iy)].re += alpha.re * x[@intCast(ix)].re - alpha.im * x[@intCast(ix)].im;
                    y[@intCast(iy)].im += alpha.re * x[@intCast(ix)].im + alpha.im * x[@intCast(ix)].re;
                    y[@intCast(iy + incy)].re += alpha.re * x[@intCast(ix + incx)].re - alpha.im * x[@intCast(ix + incx)].im;
                    y[@intCast(iy + incy)].im += alpha.re * x[@intCast(ix + incx)].im + alpha.im * x[@intCast(ix + incx)].re;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)].re += alpha.re * x[@intCast(ix)].re - alpha.im * x[@intCast(ix)].im;
                y[@intCast(iy)].im += alpha.re * x[@intCast(ix)].im + alpha.im * x[@intCast(ix)].re;

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.axpy only supports simple types."),
        .unsupported => unreachable,
    }
}
