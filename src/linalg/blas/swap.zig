const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const options = @import("options");

pub inline fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    switch (numericType) {
        .bool => @compileError("blas.swap does not support bool."),
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
                    const x0 = x[@intCast(ix)];
                    const x1 = x[@intCast(ix + incx)];
                    const x2 = x[@intCast(ix + incx2)];
                    const x3 = x[@intCast(ix + incx3)];

                    x[@intCast(ix)] = y[@intCast(iy)];
                    y[@intCast(iy)] = x0;
                    x[@intCast(ix + incx)] = y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = x1;
                    x[@intCast(ix + incx2)] = y[@intCast(iy + incy2)];
                    y[@intCast(iy + incy2)] = x2;
                    x[@intCast(ix + incx3)] = y[@intCast(iy + incy3)];
                    y[@intCast(iy + incy3)] = x3;

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const x0 = x[@intCast(ix)];

                x[@intCast(ix)] = y[@intCast(iy)];
                y[@intCast(iy)] = x0;

                ix += incx;
                iy += incy;
            }
        },
        .cfloat => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    const x0 = x[@intCast(ix)];
                    const x1 = x[@intCast(ix + incx)];

                    x[@intCast(ix)] = y[@intCast(iy)];
                    y[@intCast(iy)] = x0;
                    x[@intCast(ix + incx)] = y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = x1;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const x0 = x[@intCast(ix)];

                x[@intCast(ix)] = y[@intCast(iy)];
                y[@intCast(iy)] = x0;

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.swap only supports simple types."),
        .unsupported => unreachable,
    }
}
