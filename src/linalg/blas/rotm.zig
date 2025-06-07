const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: [*]const T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n == 0 or param[0] == -2) return;

    switch (numericType) {
        .bool => @compileError("blas.rotm does not support bool."),
        .int, .float => {
            const flag = param[0];

            if (flag == -1) {
                const h11 = param[1];
                const h21 = param[2];
                const h12 = param[3];
                const h22 = param[4];
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
                        x[@intCast(ix)] = h11 * x[@intCast(ix)] + h12 * y[@intCast(iy)];
                        y[@intCast(iy)] = h21 * x0 + h22 * y[@intCast(iy)];
                        x[@intCast(ix + incx)] = h11 * x[@intCast(ix + incx)] + h12 * y[@intCast(iy + incy)];
                        y[@intCast(iy + incy)] = h21 * x1 + h22 * y[@intCast(iy + incy)];
                        x[@intCast(ix + incx2)] = h11 * x[@intCast(ix + incx2)] + h12 * y[@intCast(iy + incy2)];
                        y[@intCast(iy + incy2)] = h21 * x2 + h22 * y[@intCast(iy + incy2)];
                        x[@intCast(ix + incx3)] = h11 * x[@intCast(ix + incx3)] + h12 * y[@intCast(iy + incy3)];
                        y[@intCast(iy + incy3)] = h21 * x3 + h22 * y[@intCast(iy + incy3)];

                        ix += incx4;
                        iy += incy4;
                    }
                }

                for (@intCast(nu)..@intCast(n)) |_| {
                    const x0 = x[@intCast(ix)];
                    x[@intCast(ix)] = h11 * x[@intCast(ix)] + h12 * y[@intCast(iy)];
                    y[@intCast(iy)] = h21 * x0 + h22 * y[@intCast(iy)];

                    ix += incx;
                    iy += incy;
                }
            } else if (flag == 0) {
                const h21 = param[2];
                const h12 = param[3];
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
                        x[@intCast(ix)] = x[@intCast(ix)] + h12 * y[@intCast(iy)];
                        y[@intCast(iy)] = y[@intCast(iy)] + h21 * x0;
                        x[@intCast(ix + incx)] = x[@intCast(ix + incx)] + h12 * y[@intCast(iy + incy)];
                        y[@intCast(iy + incy)] = y[@intCast(iy + incy)] + h21 * x1;
                        x[@intCast(ix + incx2)] = x[@intCast(ix + incx2)] + h12 * y[@intCast(iy + incy2)];
                        y[@intCast(iy + incy2)] = y[@intCast(iy + incy2)] + h21 * x2;
                        x[@intCast(ix + incx3)] = x[@intCast(ix + incx3)] + h12 * y[@intCast(iy + incy3)];
                        y[@intCast(iy + incy3)] = y[@intCast(iy + incy3)] + h21 * x3;

                        ix += incx4;
                        iy += incy4;
                    }
                }

                for (@intCast(nu)..@intCast(n)) |_| {
                    const x0 = x[@intCast(ix)];
                    x[@intCast(ix)] = x[@intCast(ix)] + h12 * y[@intCast(iy)];
                    y[@intCast(iy)] = y[@intCast(iy)] + h21 * x0;

                    ix += incx;
                    iy += incy;
                }
            } else if (flag == 1) {
                const h11 = param[1];
                const h22 = param[4];
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
                        x[@intCast(ix)] = h11 * x[@intCast(ix)] + y[@intCast(iy)];
                        y[@intCast(iy)] = h22 * y[@intCast(iy)] - x0;
                        x[@intCast(ix + incx)] = h11 * x[@intCast(ix + incx)] + y[@intCast(iy + incy)];
                        y[@intCast(iy + incy)] = h22 * y[@intCast(iy + incy)] - x1;
                        x[@intCast(ix + incx2)] = h11 * x[@intCast(ix + incx2)] + y[@intCast(iy + incy2)];
                        y[@intCast(iy + incy2)] = h22 * y[@intCast(iy + incy2)] - x2;
                        x[@intCast(ix + incx3)] = h11 * x[@intCast(ix + incx3)] + y[@intCast(iy + incy3)];
                        y[@intCast(iy + incy3)] = h22 * y[@intCast(iy + incy3)] - x3;

                        ix += incx4;
                        iy += incy4;
                    }
                }

                for (@intCast(nu)..@intCast(n)) |_| {
                    const x0 = x[@intCast(ix)];
                    x[@intCast(ix)] = h11 * x[@intCast(ix)] + y[@intCast(iy)];
                    y[@intCast(iy)] = h22 * y[@intCast(iy)] - x0;

                    ix += incx;
                    iy += incy;
                }
            }
        },
        .cfloat => @compileError("blas.rotm does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.rotm only supports simple types."),
    }
}
