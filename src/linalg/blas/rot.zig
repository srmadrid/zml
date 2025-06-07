const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

const Scalar = types.Scalar;

pub inline fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Scalar(T), s: Scalar(T)) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or (c == 1 and s == 0)) return;

    switch (numericType) {
        .bool => @compileError("blas.rot does not support bool."),
        .int, .float => {
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
                    var temp = c * x[@intCast(ix)] + s * y[@intCast(iy)];
                    y[@intCast(iy)] = c * y[@intCast(iy)] - s * x[@intCast(ix)];
                    x[@intCast(ix)] = temp;

                    temp = c * x[@intCast(ix + incx)] + s * y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = c * y[@intCast(iy + incy)] - s * x[@intCast(ix + incx)];
                    x[@intCast(ix + incx)] = temp;

                    temp = c * x[@intCast(ix + incx2)] + s * y[@intCast(iy + incy2)];
                    y[@intCast(iy + incy2)] = c * y[@intCast(iy + incy2)] - s * x[@intCast(ix + incx2)];
                    x[@intCast(ix + incx2)] = temp;

                    temp = c * x[@intCast(ix + incx3)] + s * y[@intCast(iy + incy3)];
                    y[@intCast(iy + incy3)] = c * y[@intCast(iy + incy3)] - s * x[@intCast(ix + incx3)];
                    x[@intCast(ix + incx3)] = temp;

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const temp = c * x[@intCast(ix)] + s * y[@intCast(iy)];
                y[@intCast(iy)] = c * y[@intCast(iy)] - s * x[@intCast(ix)];
                x[@intCast(ix)] = temp;

                ix += incx;
                iy += incy;
            }
        },
        .cfloat => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    var temp = c * x[@intCast(ix)].re + s * y[@intCast(iy)].re;
                    y[@intCast(iy)].re = c * y[@intCast(iy)].re - s * x[@intCast(ix)].re;
                    x[@intCast(ix)].re = temp;

                    temp = c * x[@intCast(ix)].im + s * y[@intCast(iy)].im;
                    y[@intCast(iy)].im = c * y[@intCast(iy)].im - s * x[@intCast(ix)].im;
                    x[@intCast(ix)].im = temp;

                    temp = c * x[@intCast(ix + incx)].re + s * y[@intCast(iy + incy)].re;
                    y[@intCast(iy + incy)].re = c * y[@intCast(iy + incy)].re - s * x[@intCast(ix + incx)].re;
                    x[@intCast(ix + incx)].re = temp;

                    temp = c * x[@intCast(ix + incx)].im + s * y[@intCast(iy + incy)].im;
                    y[@intCast(iy + incy)].im = c * y[@intCast(iy + incy)].im - s * x[@intCast(ix + incx)].im;
                    x[@intCast(ix + incx)].im = temp;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const temp = c * x[@intCast(ix)].re + s * y[@intCast(iy)].re;
                y[@intCast(iy)].re = c * y[@intCast(iy)].re - s * x[@intCast(ix)].re;
                x[@intCast(ix)].re = temp;

                const temp2 = c * x[@intCast(ix)].im + s * y[@intCast(iy)].im;
                y[@intCast(iy)].im = c * y[@intCast(iy)].im - s * x[@intCast(ix)].im;
                x[@intCast(ix)].im = temp2;

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.rot only supports simple types."),
    }
}
